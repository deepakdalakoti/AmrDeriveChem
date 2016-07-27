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
#include "ChemDriver.H"
#include "WritePlotFile.H"

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
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

    if (pp.contains("verbose"))
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

    ChemDriver cd;

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idYin = -1;
    int idTin = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "Y(" + cd.speciesNames()[0] + ")";
    const std::string TName = "temp";
    const int nSpec = cd.numSpecies();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == spName) idYin = i;
        if (plotVarNames[i] == TName) idTin = i;
    }
    if (ParallelDescriptor::IOProcessor() && (idYin<0 || idTin<0)) {
        BoxLib::Abort("Cannot find required data in pltfile");
    }

    const int idCout = 0;
    const int idTout = nSpec;
    const int nCompIn = nSpec+1;
    const int nCompOut = nSpec + 1;
    Real Patm = 1.0; pp.query ("Patm", Patm);
    
    Array<std::string> outNames(nCompOut);
    Array<std::string> inNames(nCompIn);
    Array<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Ys start here
    const int idTlocal = nSpec; // T starts here
    for (int i=0; i<nSpec; ++i)
    {
        destFillComps[i] = idYlocal + i;
        inNames[i] =  "Y(" + cd.speciesNames()[i] + ")";
        outNames[i] = "C(" + cd.speciesNames()[i] + ")";
    }
    destFillComps[idTlocal] = idTlocal;
    inNames[idTlocal] = TName;
    outNames[idTout] = TName;
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "To read:" << std::endl;
        for (int i=0; i<nCompIn; ++i) {
	  std::cout << inNames[i] << " into comp: " << destFillComps[i] << std::endl;
	}
        std::cout << "To write:" << std::endl;
        for (int i=0; i<nCompOut; ++i) {
	  std::cout << outNames[i] << " into comp: " << destFillComps[i] << std::endl;
	}
	std::cout << "T from " << idTin << " in pltfile to " << idTlocal << " local, to " << idTout << " in output" << endl;
	std::cout << "Y from " << idYin << " in pltfile to " << idYlocal << " local, to " << idCout << " in output" << endl;
    }
    
    PArray<MultiFab> outdata(Nlev,PArrayManage);
    const int nGrow = 0;

    FArrayBox Yt;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        outdata.set(lev,new MultiFab(ba,nCompOut,nGrow));
        MultiFab indata(ba,nCompIn,nGrow);

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        amrData.FillVar(indata,lev,inNames,destFillComps);
        for (int i=0; i<inNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        for (MFIter mfi(indata); mfi.isValid(); ++mfi)
        {
            const FArrayBox& Ylocal = indata[mfi];
            const FArrayBox& Tlocal = indata[mfi];
            FArrayBox& Cout = outdata[lev][mfi];
            FArrayBox& Tout = outdata[lev][mfi];
            const Box& box = mfi.validbox();

            Yt.resize(box, nSpec);
	    int idYtmp = 0;
            cd.normalizeMassFrac(Yt,Ylocal,"N2",box,idYlocal,0);
            cd.massFracToMolarConc (Cout, Yt, Tlocal, Patm, box, idYtmp, idTlocal, idCout);
            Tout.copy(Tlocal,idTlocal,idTout,1);
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive finished for level " << lev << endl;
    }

    std::string nfile(getFileRoot(plotFileName) + "_C");

    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << nfile << endl;
    
    const bool verb = false;
    const AmrData& a = amrData;
    WritePlotfile("NavierStokes-V1.1",outdata,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),nfile,outNames,verb);

    BoxLib::Finalize();
    return 0;
}
