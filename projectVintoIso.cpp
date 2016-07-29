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
#include "projectVintoIso_F.H"
#include "WritePlotFile.H"
#include "AppendToPlotFile.H"

static void print_usage (int argc,
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

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool verbose=false;
    if (pp.contains("verbose"))
    {
        verbose=true;
        AmrData::SetVerbose(true);
    }

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    std::string progressName = "temp"; pp.query("progressName",progressName);
    Array<int> vComps(BL_SPACEDIM,-1);
    int nc = pp.countval("vel");
    if (nc!=BL_SPACEDIM)
        BoxLib::Abort("must specify a velocity component for each spatial dimension!");
    pp.getarr("vel",vComps,0,nc);

    int idC = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == progressName) idC = i;
    }
    if (ParallelDescriptor::IOProcessor())
    {
        if (idC<0)
            cerr << "Cannot find required data in pltfile" << endl;
        for (int i=0; i<BL_SPACEDIM; ++i)
            if (vComps[i]<0)
                cerr << "Cannot find vector component in pltfile" << endl;
    }

    const int idCst = 0;
    const int idVst = idCst + 1;
    const int nCompIn = idVst + BL_SPACEDIM;

    Array<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];
    for (int i=0; i<BL_SPACEDIM; ++i)
        inVarNames[idVst+i] = plotVarNames[vComps[i]];
    if (ParallelDescriptor::IOProcessor() && verbose)
        for (int i=0; i<inVarNames.size(); ++i)
            std::cout << "going to read: " << inVarNames[i] << endl;

    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;

    const int idPV = nCompIn;
    const int nCompOut = idPV + BL_SPACEDIM;

    PArray<MultiFab> state(Nlev,PArrayManage);
    PArray<Geometry> geoms(Nlev,PArrayManage);
    const int nGrow = 1;

    FArrayBox tmp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        state.set(lev,new MultiFab(ba,nCompOut,nGrow));
        const Array<Real>& delta = amrData.DxLevel()[lev];

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        state[lev].setBndry(0.0); // Initialize grow cells to 0....here don't care to get grads right on c-f
        amrData.FillVar(state[lev],lev,inVarNames,destFillComps);
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        const bool do_corners = false;
        geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));

        // Fix up grow cells.  Use extrap for guess
        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                          s.dataPtr(idCst),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                          &nCompIn);
        }

        // Fix up fine-fine and periodic for idCst
        state[lev].FillBoundary(idCst,nCompIn);
        geoms[lev].FillPeriodicBoundary(state[lev],idCst,nCompIn,do_corners);

        // Project vector into progress contours
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_PRJECTV(box.loVect(),box.hiVect(),
                         s.dataPtr(idCst), ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         s.dataPtr(idVst), ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         s.dataPtr(idPV),  ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         delta.dataPtr());
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Vector has been projected on level " << lev << endl;

    }

    Array<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i)
        nnames[i] = inVarNames[i];
    nnames[idPV+0] = plotVarNames[vComps[0]] + "_pgx";
    nnames[idPV+1] = plotVarNames[vComps[1]] + "_pgy";
#if BL_SPACEDIM==3
    nnames[idPV+2] = plotVarNames[vComps[2]] + "_pgz";
#endif

    bool verb=false;
    bool appendPlotFile=false; pp.query("appendPlotFile",appendPlotFile);
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
        std::string newHeaderName = "NEWHEADER"; pp.query("newHeaderName",newHeaderName);
        AppendToPlotFile(amrData,ostate,plotFileName,namesOut,newMFBaseName,newHeaderName,verb);
        
        if (ParallelDescriptor::IOProcessor())
        {
            cerr << "...finished.  Note: to see new data, you must rename NewHeader in the" << endl;
            cerr << "              pltfile to Header (probably want to save the original first)" << endl;
        }
    }
    else
    {
        std::string nfile(getFileRoot(plotFileName) + "_proj");

        if (ParallelDescriptor::IOProcessor())
            cout << "Writing new data to " << nfile << endl;

        const AmrData& a = amrData;
        WritePlotfile("NavierStokes-V1.1",state,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                      a.DxLevel(),a.CoordSys(),nfile,nnames,verb);
    }

    BoxLib::Finalize();
    return 0;
}
