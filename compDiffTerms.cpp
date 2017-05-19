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
#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include "compDiffTerms_F.H"

#include "WritePlotFile.H"

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
    std::cerr << "\t TransportFile=[string]\n";
    std::cerr << "\t finestLevel=[int]\n";
    std::cerr << "\t verbose=[int]\n";
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

   ChemDriver cd;

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
   Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();


    int idXin = -1;
    int idTin = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "Y(" + cd.speciesNames()[0] + ")";
    const int nSpec = cd.numSpecies();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == spName) idXin = i;
        if (plotVarNames[i] == "temp") idTin = i;
    }
    if (ParallelDescriptor::IOProcessor() && (idXin<0 || idTin<0) )
        cerr << "Cannot find required data in pltfile" << endl;


    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    const int idXst = 0;
    const int idTst = nSpec;
    const int nCompIn = idTst + 1;

    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;
    Array<std::string> inVarNames(nCompIn);
    for (int i=0; i<nSpec; ++i)
        inVarNames[idXst+i] = plotVarNames[idXin+i];
    inVarNames[idTst] = plotVarNames[idTin];

    PArray<MultiFab> outState(Nlev,PArrayManage);
    PArray<Geometry> geoms(Nlev,PArrayManage);
    const int nGrow = 1;

    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        MultiFab inState(ba,nCompIn,nGrow);
        outState.set(lev,new MultiFab(ba,nSpec+1,0));

        Array<Real> dx(BL_SPACEDIM);
        for (int i=0; i<BL_SPACEDIM; ++i)
            dx[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        amrData.FillVar(inState,lev,inVarNames,destFillComps);
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        const bool do_corners = false;
        geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));

        // Extrap grow cells, as a first guess
        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(inState); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = inState[mfi];
            const Box& box = mfi.validbox();
            FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                          s.dataPtr(),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                          &nCompIn);
        }

        // Fix up fine-fine and periodic
        inState.FillBoundary(0,nCompIn);
        geoms[lev].FillPeriodicBoundary(inState,0,nCompIn,do_corners);

        FArrayBox Y,H,RhoD;
        for (MFIter mfi(inState); mfi.isValid(); ++mfi)
        {
            const Box& vbox = mfi.validbox();
            const Box gbox = Box(vbox).grow(nGrow);
            H.resize(gbox,nSpec);
            Y.resize(gbox,nSpec);
            RhoD.resize(gbox,nSpec+1); // last comp will hold lambda
            const FArrayBox& X = inState[mfi];
            const FArrayBox& T = inState[mfi];

            FArrayBox& DTerms = outState[lev][mfi];

            // Compute some stuff on grown box
            cd.moleFracToMassFrac(Y,X,gbox,idXst,0);
            cd.getHGivenT(H,T,gbox,idTst,0);

            const Real Patm = 60.0; 
            const int do_temp = 1;
            const int do_VelVisc = 0;

#if 1
            FORT_MIXAVG_RHODIFF_TEMP(gbox.loVect(),gbox.hiVect(),
                                     RhoD.dataPtr(),ARLIM(RhoD.loVect()),ARLIM(RhoD.hiVect()),
                                     T.dataPtr(idTst),ARLIM(T.loVect()),ARLIM(T.hiVect()),
                                     Y.dataPtr(),ARLIM(Y.loVect()),ARLIM(Y.hiVect()),
                                     &Patm, &do_temp, &do_VelVisc);
#else
            RhoD.setVal(0.0);
#endif
#if 1
            FORT_COMPDIFFTERMS(vbox.loVect(),vbox.hiVect(),
                               T.dataPtr(idTst),    ARLIM(T.loVect()),      ARLIM(T.hiVect()),
                               Y.dataPtr(),         ARLIM(Y.loVect()),      ARLIM(Y.hiVect()),
                               H.dataPtr(),         ARLIM(H.loVect()),      ARLIM(H.hiVect()),
                               RhoD.dataPtr(),      ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               RhoD.dataPtr(nSpec), ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               DTerms.dataPtr(),    ARLIM(DTerms.loVect()) ,ARLIM(DTerms.hiVect()),
                               dx.dataPtr(), &nSpec);
#else
            DTerms.setVal(0.0);
#endif
        }        

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive completed on level " << lev << endl;

    }

    std::string outfile(getFileRoot(plotFileName) + "_diff"); pp.query("outfile",outfile);
    
    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << outfile << endl;
    
    Array<string> nnames(nSpec+1);
    nnames[nSpec] = "Div(lambda.Grad(T)) + Div(rho.D_i.h_i.Grad(Y_i))";
    for (int i=0; i<nSpec; ++i) {
        nnames[i] = "Div(rho.D_" + cd.speciesNames()[i] + ".Grad(Y_" + cd.speciesNames()[i] + "))";
    }
    bool verb=false;
    const AmrData& a = amrData;
    WritePlotfile("NavierStokes-V1.1",outState,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),outfile,nnames,verb);

    BoxLib::Finalize();
    return 0;
}
