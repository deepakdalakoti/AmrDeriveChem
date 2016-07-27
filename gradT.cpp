//BL_COPYRIGHT_NOTICE
#include <winstd.H>

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>

using std::cout;
using std::cerr;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "FArrayBox.H"
#include "Utility.H"

#include "Geometry.H"
#include "gradT_F.H"

#include <AppendToPlotFile.H>
#include <WritePlotFile.H>

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile plotFileNames=f1 f2 <where time(f1)<time(f2)> [options] \n\tOptions:\n";
    exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
    vector<std::string> tokens = BoxLib::Tokenize(infile,std::string("/"));
    return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    const Real strt_time = ParallelDescriptor::second();
    Real strt_io, io_time = 0;

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    int verbose=0; pp.query("verbose",verbose);
    if (verbose>1)
        AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    strt_io = ParallelDescriptor::second();
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();
    io_time += ParallelDescriptor::second() - strt_io;

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    std::string progressName = "temp"; pp.query("progressName",progressName);

    int idC = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == progressName) idC = i;
    }
    if (ParallelDescriptor::IOProcessor() && idC<0) {
        BoxLib::Abort("Cannot find required data in pltfile");
    }

    const int idCst = 0;
    const int nCompIn = idCst + 1;

    Array<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];
    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;
    

    const int idGr = nCompIn;
    const int nCompOut = idGr + BL_SPACEDIM +1 ; // 1 component stores the ||gradT||

    PArray<MultiFab> state(Nlev,PArrayManage);
    PArray<Geometry> geoms(Nlev,PArrayManage);
    const int nGrow = 1;

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "pltfile has Nlev: " << Nlev << std::endl;
    }
    for (int lev=0; lev<Nlev; ++lev)
    {
        strt_io = ParallelDescriptor::second();
        const BoxArray ba = amrData.boxArray(lev);
        state.set(lev,new MultiFab(ba,nCompOut,nGrow));
        const Array<Real>& delta = amrData.DxLevel()[lev];

        // Get input state data onto intersection ba
        const int myNComp = 1; // gonna need this for fortran calls

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        state[lev].setBndry(0.0); // Initialize grow cells to 0....here don't care to get grads right on c-f

        ParallelDescriptor::Barrier();
#if 0
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_HACKVAL(box.loVect(), box.hiVect(),
                         s.dataPtr(idCst),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         delta.dataPtr(), amrData.ProbLo().dataPtr());
        }
#else
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Doing FillVar at lev: " << lev << std::endl;
        }
        amrData.FillVar(state[lev],lev,inVarNames,destFillComps);

        //for (int i=0; i<inVarNames.size(); ++i)
        //  amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));
#endif
        ParallelDescriptor::Barrier();

        io_time += ParallelDescriptor::second() - strt_io;
        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        const bool do_corners = false;
        geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));

#if 1
        // Fix up grow cells.  Use extrap for guess
        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                          s.dataPtr(idCst),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                          &myNComp);
        }
#endif

        // Fix up fine-fine and periodic for idCst
        state[lev].FillBoundary(idCst,1);
        geoms[lev].FillPeriodicBoundary(state[lev],idCst,1,do_corners);

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Boundary cells filled at lev: " << lev << std::endl;
        }
        // Compute mag of progress gradient.  Result in state, comp=idGr
        FArrayBox nWork;
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_GRAD(box.loVect(),box.hiVect(),
                      s.dataPtr(idCst),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                      s.dataPtr(idGr), ARLIM(s.loVect()),ARLIM(s.hiVect()),
                      delta.dataPtr());
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Gradient vector has been computed on level " << lev << endl;

        ParallelDescriptor::Barrier();
    }

    /* adding one component in the nCompOut to get the ||grad progress ||, and this might be useful to compute the thermal thickness. */
    Array<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i)
       nnames[i] = inVarNames[i];
    nnames[idGr+0] = progressName + "_gx";
    nnames[idGr+1] = progressName + "_gy";
#if BL_SPACEDIM==3
    nnames[idGr+2] = progressName + "_gz";
#endif
    nnames[idGr+BL_SPACEDIM] = "||grad"+ progressName+ "||";
    bool verb=false;
    bool appendPlotFile=false; pp.query("appendPlotFile",appendPlotFile);

    strt_io = ParallelDescriptor::second();
    if (appendPlotFile)
       {

        int nStateOut = nCompOut - nCompIn;
        PArray<MultiFab> ostate(Nlev,PArrayManage);
        for (int lev=0; lev<Nlev; ++lev)
        {
            const BoxArray ba = state[lev].boxArray();
            ostate.set(lev,new MultiFab(ba,nStateOut,nGrow));
            MultiFab::Copy(ostate[lev],state[lev],nCompIn,0,nStateOut,0);
        }
        Array<std::string> namesOut(nStateOut);
        for (int i=0; i<nStateOut; ++i)
            namesOut[i] = nnames[nCompIn+i];

        std::string newMFBaseName = "NEWDAT"; pp.query("newMFBaseName",newMFBaseName);
        std::string newHeaderName = "NewHeader"; pp.query("newHeaderName",newHeaderName);
        ParallelDescriptor::Barrier();
        AppendToPlotFile(amrData,ostate,plotFileName,namesOut,newMFBaseName,newHeaderName,verb);
        ParallelDescriptor::Barrier();
        
        if (ParallelDescriptor::IOProcessor())
        {
            cerr << "...finished.  Note: to see new data, you must rename NewHeader in the" << endl;
            cerr << "              pltfile to Header (probably want to save the original first)" << endl;
        }
    }
    else
    {
        std::string outfile(getFileRoot(plotFileName) + "_gt"); pp.query("outfile",outfile);

        if (ParallelDescriptor::IOProcessor())
            cout << "Writing new data to " << outfile << endl;

        const AmrData& a = amrData;
        ParallelDescriptor::Barrier();
        WritePlotfile("NavierStokes-V1.1",state,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                      a.DxLevel(),a.CoordSys(),outfile,nnames,verb);
        ParallelDescriptor::Barrier();
    }

    io_time += ParallelDescriptor::second() - strt_io;

    //long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
    //long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
    
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
/*
    ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nFAB byte spread across MPI nodes: ["
                  << min_fab_bytes
                  << " ... "
                  << max_fab_bytes
                  << "]\n";
    }
*/
    const Real end_time = ParallelDescriptor::second();
    const Real run_time = end_time - strt_time - io_time;

    Real run_time_max, run_time_min; run_time_max = run_time_min = run_time;
    Real io_time_max, io_time_min; io_time_max = io_time_min = io_time;

    ParallelDescriptor::ReduceRealMax(run_time_max,IOProc);
    ParallelDescriptor::ReduceRealMax(io_time_max,IOProc);
    ParallelDescriptor::ReduceRealMin(run_time_min,IOProc);
    ParallelDescriptor::ReduceRealMin(io_time_min,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Max Compute time: " << run_time_max << '\n';
        std::cout << "Min Compute time: " << run_time_min << '\n';
        std::cout << "Max I/O time: " << io_time_max << '\n';
        std::cout << "Min I/O time: " << io_time_min << '\n';
    }


    BoxLib::Finalize();
    return 0;
}
