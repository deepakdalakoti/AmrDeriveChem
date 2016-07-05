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
#include "ArrayLim.H"
#include "summassCyl_F.H"

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=f1 [options] \n\tOptions:\n";
    exit(1);
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

    int verbose=0;
    pp.query("verbose",verbose);
    if (verbose>1)
        AmrData::SetVerbose(true);

    std::string infile; pp.get("infile",infile);
    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);

    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    std::string compName; pp.get("compName",compName);

    Real rMax=0.; pp.query("rMax",rMax);
    Real zMax=0.; pp.query("zMax",zMax);

    Real mySum = 0;
    const Array<Real>& plo = amrData.ProbLo();
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray& ba = amrData.boxArray(lev);
        MultiFab data(ba,1,0,Fab_allocate);
        amrData.FillVar(data,lev,compName);
        if (ParallelDescriptor::IOProcessor())
            cout << "...data filled on level " << lev << '\n';

        // Zero the covered data
        if (lev<finestLevel)
        {
	    BoxArray baf = BoxArray(amrData.boxArray(lev+1)).coarsen(amrData.RefRatio()[lev]);	  

            for (MFIter mfi(data); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.validbox();
                FArrayBox& fab = data[mfi];
                std::vector< std::pair<int,Box> > isects = baf.intersections(box);
            
                for (int i = 0; i < isects.size(); i++)
                    fab.setVal(0,isects[i].second,0,1);
            }
        }

        // volume-weighted sum
        Real cellVol = 1;
        Array<Real> dx(BL_SPACEDIM);
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            dx[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);
            cellVol *= dx[i];
        }
        data.mult(cellVol);

        for (MFIter mfi(data); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const FArrayBox& fab = data[mfi];
            mySum += FORT_SUM_IN_RZ(box.loVect(), box.hiVect(),
                                    fab.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()),
                                    dx.dataPtr(), plo.dataPtr(), &rMax, &zMax);
        }
    }

    ParallelDescriptor::ReduceRealSum(mySum);

    if (ParallelDescriptor::IOProcessor())
        cout << amrData.Time() << " " << mySum << '\n';

    BoxLib::Finalize();
    return 0;
}
