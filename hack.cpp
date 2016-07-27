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
#include "hack_F.H"
#include "WritePlotFile.H"
#include "AppendToPlotFile.H"

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

    int idC = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == progressName) idC = i;
    }
    if (ParallelDescriptor::IOProcessor() && idC<0)
        cerr << "Cannot find required data in pltfile" << endl;

    const int idCst = 0;
    const int nCompIn = idCst + 1;

    Array<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];
    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;
    
    const int nCompOut = nCompIn;

    PArray<MultiFab> state(Nlev,PArrayManage);
    const int nGrow = 0;

    FArrayBox tmp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        state.set(lev,new MultiFab(ba,nCompOut,nGrow));
        Array<Real> delta(BL_SPACEDIM);
        for (int i=0; i<BL_SPACEDIM; ++i)
            delta[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);


        // Get input state data onto intersection ba
        const int myNComp = 1; // gonna need this for fortran calls

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        if (nGrow>0)
            state[lev].setBndry(0.0); // Initialize grow cells to 0....here don't care to get grads right on c-f
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_HACKVAL(box.loVect(), box.hiVect(),
                         s.dataPtr(idCst),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         delta.dataPtr(), amrData.ProbLo().dataPtr());
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been filled for level " << lev << endl;

    }

    Array<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i)
        nnames[i] = inVarNames[i];

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
        std::string outfile(getFileRoot(plotFileName) + "_h"); pp.query("outfile",outfile);

        if (ParallelDescriptor::IOProcessor())
            cout << "Writing new data to " << outfile << endl;

        const AmrData& a = amrData;
        WritePlotfile("NavierStokes-V1.1",state,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                      a.DxLevel(),a.CoordSys(),outfile,nnames,verb);
    }

    BoxLib::Finalize();
    return 0;
}
