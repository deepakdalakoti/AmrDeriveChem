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

typedef ChemDriver::Edge Edge;
typedef std::list<Edge> EdgeList;

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
    int nReac = cd.numReactions();
    int sCompT = nReac;
    int sCompR = nReac+1;
    int sCompY = 0;
    Array<string> varNames(nReac+2);
    for (int i=0; i<nReac; ++i) {
      varNames[sCompY + i] = BoxLib::Concatenate("R",i+1,5);
    }
    varNames[sCompT] = "temp";
    varNames[sCompR] = "density";

    Array<int> destFillComps(varNames.size());
    for (int i=0; i<varNames.size(); ++i)
        destFillComps[i] = i;
    
    Array<std::string> species;
    pp.getarr("species",species,0,pp.countval("species"));

    std::string QPDatom="C"; pp.query("QPDatom",QPDatom);

    EdgeList edges = cd.getEdges(QPDatom);
    std::vector<std::map<int,Real> > rlist;
    for (int j=0; j<species.size(); ++j) {
      for (EdgeList::const_iterator it=edges.begin(); it!=edges.end(); ++it) {
	if (it->touchesSp(species[j])) {
	  const std::string& spL = it->left();
	  const std::string& spR = it->right();
	  const Array<std::pair<int,Real> > RWL = it->RateWeightList();
          std::map<int,Real> r;
	  for (int i=0; i<RWL.size(); ++i) {
	    r[RWL[i].first] = ( spL == species[j]  ?  -1  : +1 ) * RWL[i].second;
	  }
          rlist.push_back(r);
	}
      }
    }

    PArray<MultiFab> CD(Nlev);
    FArrayBox RhoH, RhoYdot;
    for (int iLevel=0; iLevel<Nlev; ++iLevel)
    {
        CD.set(iLevel,new MultiFab(amrData.boxArray(iLevel),species.size()+1,0,Fab_allocate));
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

	FArrayBox tmp;
 static bool first = true;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
	  FArrayBox& cd = CD[iLevel][mfi];
	  tmp.resize(cd.box(),1);
	  cd.setVal(0);
	  for (int i=0; i<species.size(); ++i) {
	    for (std::map<int,Real>::const_iterator rit = rlist[i].begin(); rit!=rlist[i].end(); ++rit) {
	      int reacComp = rit->first;
	      Real reacMult = rit->second;
  if (first && ParallelDescriptor::IOProcessor()) {
        std::cout << "s,r,c: " << species[i] << " " << reacComp << " " << reacMult << std::endl;
}
	      tmp.copy(mf[mfi],reacComp,0,1);
	      tmp.mult(reacMult);
	      cd.plus(tmp,0,i,1);
	    }
	  }
         first = false;
	}
        MultiFab::Copy(CD[iLevel],mf,sCompT,species.size(),1,0);
    }
    Array<std::string> outNames(species.size()+1);
    for (int i=0; i<species.size(); ++i) {
      outNames[i] = "C(" + species[i] + ")_Dot";
    }
    outNames[species.size()] = "temp";
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

