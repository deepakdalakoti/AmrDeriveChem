
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
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
    
    ParmParse pp;

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    
    int verbose = 0; pp.query("verbose",verbose);
    if (verbose > 2)
        AmrData::SetVerbose(true);

    std::string infile; pp.get("infile",infile);
    std::string outfile = infile + std::string("_reacs"); pp.query("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);    
    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Real Patm=60; pp.query("Patm",Patm);

    ChemDriver cd;

    int nSpec = cd.numSpecies();
    int sCompX = 0;
    int sCompT = nSpec;
    Array<string> varNames(nSpec+1);
    for (int i=0; i<nSpec; ++i)
        varNames[sCompX + i] = "Y(" + cd.speciesNames()[i] + ")";
    varNames[sCompT] = "temp";

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
    
    PArray<MultiFab> reacs(Nlev), Out(Nlev);
    int Nreacs = cd.numReactions();
    FArrayBox fwd,rev;
    Array<int> reacIds(Nreacs);
    for (int i=0; i<Nreacs; ++i)
    {
        reacIds[i] = i;
    }
    for (int iLevel=0; iLevel<Nlev; ++iLevel)
    {
        reacs.set(iLevel,new MultiFab(amrData.boxArray(iLevel),Nreacs,0,Fab_allocate));
        Out.set(iLevel,new MultiFab(amrData.boxArray(iLevel),1,0,Fab_allocate));

        MultiFab mf(amrData.boxArray(iLevel),varNames.size(),0,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);
  

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;
    
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& Y = mf[mfi];
            FArrayBox X;
            Real small = 1e-5;
            const FArrayBox& T = mf[mfi];
            const Box& box = mfi.validbox();
            X.resize(box,nSpec);
            X.setVal(0.0);
            cd.massFracToMoleFrac(X,Y,box,0,0);         

            FArrayBox& r = reacs[iLevel][mfi];
            FArrayBox& out = Out[iLevel][mfi];
            fwd.resize(box,Nreacs);
            rev.resize(box,Nreacs);
            cd.fwdRevReacRatesGivenXTP(fwd,rev,reacIds,X,Y,Patm,box,sCompX,sCompT,0,0);
            out.copy(fwd,7,0,1);
            out.minus(rev,7,0,1);
            X.copy(out,0,0,1);
            X.plus(fwd,5,0,1);
            X.minus(rev,5,0,1);
            X.plus(small);
            out.divide(X,0,0,1);
            out.abs();
            // Grab auxiliary quantities
            r.copy(X,nCompOLD,Nreacs,comps.size());
        }
        if(ParallelDescriptor::IOProcessor())
             std::cout << "Derive Completed on Level " << iLevel << std::endl; 
    }

/*  Array<string> outNames(Nreacs+comps.size());
    for (int i=0; i<Nreacs; ++i)
    {        
        outNames[i] = BoxLib::Concatenate("R",i+1);;
    }
    for (int i=0; i<comps.size(); ++i)
    {
        outNames[i+Nreacs] = varNames[nCompOLD + i];
    }
*/  Array<string> outNames(1);
     outNames[0]= "Autoignition Index";
    WritePlotFile(Out,outfile,outNames,amrData);
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
