
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
    std::string outfile = infile + std::string("_dot"); pp.query("outfile",outfile);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);    
    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Array<string> varNames(2*BL_SPACEDIM);
    Array<int> vecs(2*BL_SPACEDIM);
    if (pp.countval("vecs")==vecs.size())
    {
        pp.getarr("vecs",vecs,0,vecs.size());
        for (int j=0; j<vecs.size(); ++j)
        {
            BL_ASSERT(vecs[j]>0 && vecs[j]<=amrDataNComp());
            varNames[j] = amrData.PlotVarNames()[vecs[j]];
        }
    }

    varNames.resize(varNames.size()+1);
    int wtComp; pp.get("weight",wtComp); varNames[varNames.size()-1] = amrData.PlotVarNames()[wtComp];
    Real wtThresh=1.e-8; pp.query("wtThresh",wtThresh);

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
    
    PArray<MultiFab> dot(Nlev);
    for (int iLevel=0; iLevel<Nlev; ++iLevel)
    {
        dot.set(iLevel,new MultiFab(amrData.boxArray(iLevel),1+comps.size(),0,Fab_allocate));
        MultiFab mf(amrData.boxArray(iLevel),varNames.size(),0,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;
    
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Real* res = dot[iLevel][mfi].dataPtr();
            const FArrayBox fab = mf[mfi];
            const Real* v1x = fab.dataPtr(0);
            const Real* v1y = fab.dataPtr(1);
            const Real* v1z = fab.dataPtr(2);
            const Real* v2x = fab.dataPtr(3);
            const Real* v2y = fab.dataPtr(4);
            const Real* v2z = fab.dataPtr(5);
            const Real* wt  = fab.dataPtr(6);

            const Box& box = mfi.validbox();

            const int ni = box.length(0);
            const int nj = box.length(1);
            const int nk = box.length(2);
            for (int k=0; k<nk; k++) {
                for (int j=0; j<nj; j++) {
                    for (int i=0; i<ni; i++) {
                        const int cell = (k*nj+j)*ni+i;

                        Real mag1=std::sqrt(v1x[cell]*v1x[cell] + v1y[cell]*v1y[cell] + v1z[cell]*v1z[cell]);
                        Real mag2=std::sqrt(v2x[cell]*v2x[cell] + v2y[cell]*v2y[cell] + v2z[cell]*v2z[cell]);

                        if (wt[cell]>=wtThresh)
                        {
                            res[cell] = (v1x[cell]*v2x[cell] + v1y[cell]*v2y[cell] + v1z[cell]*v2z[cell]) / std::max(epsilon, mag1*mag2);
                        }
                        else
                        {
                            res[cell] = 0;
                        }
                    } // i
                } // j
            } // k
            // Get any other data
            if (comps.size()>0)
                dot[iLevel][mfi].copy(fab,nCompOLD,1,comps.size());
        }
    }

    Array<string> outNames(1+comps.size());
    outNames[0] = "dot";
    for (int i=0; i<comps.size(); ++i)
        outNames[i+1] = varNames[nCompOLD + i];
    WritePlotFile(dot,outfile,outNames,amrData);
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
