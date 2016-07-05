
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"

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
    
    ParmParse pp;

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    
    int verbose = 0; pp.query("verbose",verbose);
    if (verbose > 2)
        AmrData::SetVerbose(true);

    std::string infile; pp.get("infile",infile);
    std::string outfile = infile + std::string("_dkappa"); pp.query("outfile",outfile);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);    
    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Array<string> varNames(3);
    pp.get("n1",varNames[0]);
    pp.get("n2",varNames[1]);
    pp.get("weight",varNames[2]);
    Real wtThresh=1.e-8; pp.query("wtThresh",wtThresh);

    int k1 = amrData.StateNumber(varNames[0]);
    int k2 = amrData.StateNumber(varNames[1]);

    Real epsilon = 1.e-12; pp.query("epsilon",epsilon);

    Array<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0; pp.query("sComp",sComp);
        int nComp = 0; pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    int nCompOLD = varNames.size();
    if (comps.size()>0)
    {
        varNames.resize( nCompOLD + comps.size() );
        for (int i=0; i<comps.size(); ++i)
        {
            varNames[nCompOLD + i] = amrData.PlotVarNames()[comps[i]];
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        for (int i=0; i<varNames.size(); ++i)
            std::cout << "Getting component: " << varNames[i] << std::endl;

    Array<int> destFillComps(varNames.size());
    for (int i=0; i<varNames.size(); ++i)
        destFillComps[i] = i;
    
    PArray<MultiFab> dKappa(Nlev);
    for (int iLevel=0; iLevel<Nlev; ++iLevel)
    {
        dKappa.set(iLevel,new MultiFab(amrData.boxArray(iLevel),1+comps.size(),0,Fab_allocate));
        MultiFab mf(amrData.boxArray(iLevel),varNames.size(),0,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;
    
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Real* res = dKappa[iLevel][mfi].dataPtr();
            const FArrayBox fab = mf[mfi];
            const Real* src1 = fab.dataPtr(0);
            const Real* src2 = fab.dataPtr(1);
            const Real*   wt = fab.dataPtr(2);
            const Box& box = mfi.validbox();

            const int ni = box.length(0);
            const int nj = box.length(1);
#if (BL_SPACEDIM==3)
            const int nk = box.length(2);
            for (int k=0; k<nk; k++) {
#endif
                for (int j=0; j<nj; j++) {
                    for (int i=0; i<ni; i++) {
#if (BL_SPACEDIM==3)
                        const int cell = (k*nj+j)*ni+i;
#else
                        const int cell = j*ni+i;
#endif
                        Real weight = wt[cell];
                        if (weight>wtThresh)
                        {
                            const Real k1 = std::abs(src1[cell]);
                            const Real k2 = std::abs(src2[cell]);
                            
                            res[cell] = 2*(k1-k2) / std::max(epsilon, k1 + k2);
                        }
                        else
                        {
                            res[cell] = 0;
                        }
                    } // i
                } // j
#if (BL_SPACEDIM==3)
            } // k
#endif
            // Get any other data
            if (comps.size()>0)
                dKappa[iLevel][mfi].copy(fab,nCompOLD,1,comps.size());
        }
    }

    Array<string> outNames(1+comps.size());
    outNames[0] = "dKappa";
    for (int i=0; i<comps.size(); ++i)
        outNames[i+1] = varNames[nCompOLD + i];
    WritePlotFile(dKappa,outfile,outNames,amrData);
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
