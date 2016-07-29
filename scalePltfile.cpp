//BL_COPYRIGHT_NOTICE
#include <winstd.H>

#include <new>
#include <iostream>

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
#include "WritePlotFile.H"

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile infile=s comps=<i1 i2 ...> vals=<r1 r2 ...>\n\tOptions: finestLevel=<i>\n";
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

    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Array<int> comps;
    Array<Real> vals;    
    int nc = pp.countval("comps");
    if (nc<1)
        print_usage(argc,argv);
    comps.resize(nc);
    pp.getarr("comps",comps,0,nc);
    pp.getarr("vals",vals,0,nc);

    int nComp = amrData.NComp();
    const Array<std::string>& names = amrData.PlotVarNames();
    Array<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i) {
        destFillComps[i] = i;
    }
    PArray<MultiFab> data(Nlev,PArrayManage);
    const int nGrow = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        data.set(lev,new MultiFab(ba,nComp,nGrow));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;
        amrData.FillVar(data[lev],lev,names,destFillComps);
        for (int i=0; i<nComp; ++i)
            amrData.FlushGrids(i);
        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        for (int j=0; j<comps.size(); ++j)
            data[lev].mult(vals[j],comps[j],1,0);

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive finished for level " << lev << endl;
    }

    std::string nfile(getFileRoot(plotFileName) + "_scale");

    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << nfile << endl;
    
    const bool verb = false;
    const AmrData& a = amrData;
    WritePlotfile("NavierStokes-V1.1",data,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),nfile,names,verb);

    BoxLib::Finalize();
    return 0;
}
