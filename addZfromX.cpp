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

    int idXin = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "X(" + cd.speciesNames()[0] + ")";
    const int nSpec = cd.numSpecies();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == spName) idXin = i;
    }
    if (ParallelDescriptor::IOProcessor() && idXin<0)
        cerr << "Cannot find required data in pltfile" << endl;

    const int nElts = cd.numElements();
    const int idZout = 0;
    const int nCompIn = nSpec;
    const int nCompOut = idZout + nElts;

    Array<std::string> outNames(nCompOut);
    Array<std::string> inNames(nCompIn);
    Array<int> destFillComps(nCompIn);
    const int idXlocal = 0; // Xs start here
    for (int i=0; i<nCompIn; ++i)
    {
        destFillComps[i] = idXlocal + i;
        inNames[i] =  "X(" + cd.speciesNames()[i] + ")";
    }
    for (int i=0; i<nElts; ++i)
    {
        outNames[i] =  "Z(" + cd.elementNames()[i] + ")";
    }
    
    PArray<MultiFab> outdata(Nlev,PArrayManage);
    const int nGrow = 0;

    Array<Array<int> > numberOfElementZinSpeciesX(nSpec);
    for (int i=0;i<nSpec; ++i)
    {
        numberOfElementZinSpeciesX[i].resize(nElts);
        for (int j=0;j<nElts; ++j)
        {
            numberOfElementZinSpeciesX[i][j]
                = cd.numberOfElementXinSpeciesY(cd.elementNames()[j],
                                                cd.speciesNames()[i]);
        }
    }

    FArrayBox tmp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        outdata.set(lev,new MultiFab(ba,nCompOut,nGrow));
        outdata[lev].setVal(0.);
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
            const FArrayBox& X = indata[mfi];
            FArrayBox& Z = outdata[lev][mfi];
            const Box& box = mfi.validbox();

            tmp.resize(box,1);

            for (int i=0;i<nSpec; ++i)
            {
                for (int j=0;j<nElts; ++j)
                {
                    tmp.copy(X,i,0,1);
                    tmp.mult(numberOfElementZinSpeciesX[i][j],0,1);
                    Z.plus(tmp,0,j,1);
                }
            }

	    tmp.setVal(0.);
	    for (int j=0;j<nElts; ++j)
	    {
		tmp.plus(Z,j,0,1);
	    }
	    tmp.invert(1);
	    for (int j=0;j<nElts; ++j)
	    {
		Z.mult(tmp,0,j,1);
	    }
	    // also Z.divide(tmp,0,j,1);
	    
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive finished for level " << lev << endl;
    }

    std::string nfile(getFileRoot(plotFileName) + "_Z");

    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << nfile << endl;
    
    const bool verb = false;
    const AmrData& a = amrData;
    WritePlotfile("NavierStokes-V1.1",outdata,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),nfile,outNames,verb);

    BoxLib::Finalize();
    return 0;
}
