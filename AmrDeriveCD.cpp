#include <iomanip>

#include "REAL.H"
#include "Box.H"
#include "PArray.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"

static
void
WritePlotFile (const PArray<MultiFab>&   mfout,
               const std::string&        ofile,
               const Array<std::string>& names,
               AmrData&                  amrData);

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    std::cout << std::setprecision(25);
    ParmParse pp;

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    
    int verbose = 0; pp.query("verbose",verbose);
    if (verbose > 2)
        AmrData::SetVerbose(true);

    std::string infile; pp.get("infile",infile);
    std::string outfile = infile + std::string("_cd"); pp.query("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);    
    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Real Patm=1; pp.query("Patm",Patm);
    ChemDriver cd;

    int nSpec = cd.numSpecies();
    int sCompT = nSpec;
    int sCompR = nSpec+1;
    int sCompY = 0;
    Array<string> varNames(nSpec+2);

    bool hasYs = amrData.CanDerive("Y(" + cd.speciesNames()[0] + ")");
    bool hasXs = amrData.CanDerive("X(" + cd.speciesNames()[0] + ")");
    BL_ASSERT(hasYs || hasXs);

    const std::string pre = hasYs ? "Y" : "X";
    for (int i=0; i<nSpec; ++i) {
        varNames[sCompY + i] = pre + "(" + cd.speciesNames()[i] + ")";
    }
    varNames[sCompT] = "temp";
    varNames[sCompR] = "density";

    Array<int> destFillComps(varNames.size());
    for (int i=0; i<varNames.size(); ++i)
        destFillComps[i] = i;
    
    Array<std::string> cdNames;
    pp.getarr("cdNames",cdNames,0,pp.countval("cdNames"));
    Array<int> cdComps(cdNames.size());
    for (int i=0; i<cdNames.size(); ++i) {
      cdComps[i] = cd.index(cdNames[i]);
      BL_ASSERT(cdComps[i] >= 0);
      if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Index of " << cdNames[i] << " is " << cdComps[i] << std::endl;
      }
    }

    PArray<MultiFab> CD(Nlev);
    FArrayBox RhoH, RhoYdot;
    for (int iLevel=0; iLevel<Nlev; ++iLevel)
    {
        CD.set(iLevel,new MultiFab(amrData.boxArray(iLevel),cdNames.size()+1,0,Fab_allocate));
        MultiFab mf(amrData.boxArray(iLevel),varNames.size(),0,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;

	long min_fab_bytes = BoxLib::TotalBytesAllocatedInFabs();
	long max_fab_bytes = BoxLib::TotalBytesAllocatedInFabsHWM();
    
	const int IOProc = ParallelDescriptor::IOProcessorNumber();
	ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
	ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
	if (ParallelDescriptor::IOProcessor()) {
	  std::cout << "\nFAB byte spread across MPI nodes: ["
		    << min_fab_bytes
		    << " ... "
		    << max_fab_bytes
		    << "]\n";
	}

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& Y = mf[mfi];

	    if (hasXs) {
	      cd.moleFracToMassFrac(Y,Y,Y.box(),sCompY,sCompY);
	    }

            const FArrayBox& T = mf[mfi];
            const FArrayBox& Rho = mf[mfi];
            const Box& box = mfi.validbox();
            FArrayBox& cdfab = CD[iLevel][mfi];

	    RhoYdot.resize(box,nSpec+1);

	    cd.reactionRateY(RhoYdot,Y,T,Patm,box,sCompY,sCompT,0);

            for (int i=0; i<nSpec; ++i) {
              RhoYdot.mult(Rho,sCompR,i,1);
            }

	    cdfab.setVal(0);
	    for (int i=0; i<cdComps.size(); i++) {
	      for (IntVect iv=box.smallEnd(), be=box.bigEnd(); iv<=be; box.next(iv)) {
		Real val = RhoYdot(iv,cdComps[i])  / cd.speciesMolecWt()[cdComps[i]]; // Cdot_i = rhoYdot_i/W_i
		cdfab(iv,i) = val;
	      }
	    }
        }
        MultiFab::Copy(CD[iLevel],mf,sCompT,cdNames.size(),1,0);
    }
    Array<std::string> outNames(cdNames.size()+1);
    for (int i=0; i<cdNames.size(); ++i) {
      outNames[i] = "C(" + cdNames[i] + ")_Dot";
    }
    outNames[cdNames.size()] = "temp";
    WritePlotFile(CD,outfile,outNames,amrData);
    BoxLib::Finalize();
    return 0;
}

static
void
WritePlotFile (const PArray<MultiFab>&   mfout,
               const std::string&        ofile,
               const Array<std::string>& names,
               AmrData&                  amrData)
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Writing plotfile " << ofile << " ...\n";

    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(ofile,0755))
            BoxLib::CreateDirectoryFailed(ofile);

    ParallelDescriptor::Barrier();

    std::string ofileHeader(ofile); ofileHeader += "/Header";
  
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream os;
  
    os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    os.open(ofileHeader.c_str(), std::ios::out | std::ios::binary);

    if (os.fail()) BoxLib::FileOpenFailed(ofileHeader);

    //const int finestLevel = amrData.FinestLevel();
    const int finestLevel = mfout.size() - 1;

    os << amrData.PlotFileVersion() << '\n';
    os << names.size()              << '\n';

    for (int i = 0; i < names.size(); i++)
        os << names[i] << '\n';

    os << BL_SPACEDIM               << '\n';
    os << amrData.Time()            << '\n';
    os << finestLevel << '\n';
    for(int i = 0; i < BL_SPACEDIM; i++) os << amrData.ProbLo()[i] << ' ';
    os << '\n';
    for(int i = 0; i < BL_SPACEDIM; i++) os << amrData.ProbHi()[i] << ' ';
    os << '\n';
    for(int i = 0; i < finestLevel; i++) os << amrData.RefRatio()[i] << ' ';
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) os << amrData.ProbDomain()[i] << ' ';
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++)
    {
        for(int k = 0; k < BL_SPACEDIM; k++)
            os << amrData.DxLevel()[i][k] << ' ';
        os << '\n';
    }

    os << amrData.CoordSys() << '\n';
    os << "0\n";

    for (int lev = 0; lev <= finestLevel; ++lev)
    {
        const int nGrids = amrData.boxArray(lev).size();
        char buf[64];
        sprintf(buf, "Level_%d", lev);
    
        if (ParallelDescriptor::IOProcessor())
        {
            os << lev << ' ' << nGrids << ' ' << amrData.Time() << '\n';
            os << 0 << '\n';
    
            for (int i = 0; i < nGrids; ++i)
            {
                for(int n = 0; n < BL_SPACEDIM; n++)
                {
                    os << amrData.GridLocLo()[lev][i][n]
                       << ' '
                       << amrData.GridLocHi()[lev][i][n]
                       << '\n';
                }
            }

            std::string Level(ofile);
            Level += '/';
            Level += buf;
    
            if (!BoxLib::UtilCreateDirectory(Level, 0755))
                BoxLib::CreateDirectoryFailed(Level);
        }

        ParallelDescriptor::Barrier();

        static const std::string MultiFabBaseName("/MultiFab");
    
        std::string PathName(ofile);
        PathName += '/';
        PathName += buf;
        PathName += MultiFabBaseName;
    
        if (ParallelDescriptor::IOProcessor())
        {
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }

        BL_ASSERT(mfout[lev].nComp() == names.size());

        VisMF::Write(mfout[lev], PathName, VisMF::OneFilePerCPU);
    }

    os.close();

}
