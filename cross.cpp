
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "Geometry.H"

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
    std::string outfile = infile + std::string("_cross"); pp.query("outfile",outfile);

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
    
    const int nGrow = 1;
    const int nCompNEW = 9;
    PArray<MultiFab> cross(Nlev);
    PArray<Geometry> geoms(Nlev,PArrayManage);
    for (int iLevel=0; iLevel<Nlev; ++iLevel)
    {
        geoms.set(iLevel,new Geometry(amrData.ProbDomain()[iLevel]));
        cross.set(iLevel,new MultiFab(amrData.boxArray(iLevel),nCompNEW+comps.size(),nGrow,Fab_allocate));
        cross[iLevel].setVal(0);

        MultiFab mf(amrData.boxArray(iLevel),varNames.size(),nGrow,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;
    
        // Get auxiliary data
        if (comps.size()>0)
            MultiFab::Copy(cross[iLevel],mf,nCompOLD,nCompNEW,comps.size(),0);

        // Fix up fine-fine and periodic cells for source data (dont bother with auxiliary data)
        mf.FillBoundary(0,nCompOLD);
        geoms[iLevel].FillPeriodicBoundary(mf,0,nCompOLD);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& ofab = cross[iLevel][mfi];
            Real* outx   = ofab.dataPtr(0);
            Real* outy   = ofab.dataPtr(1);
            Real* outz   = ofab.dataPtr(2);
            Real* outmag = ofab.dataPtr(3);
            Real* fx     = ofab.dataPtr(4);
            Real* fy     = ofab.dataPtr(5);
            Real* fz     = ofab.dataPtr(6);
            Real* fmag   = ofab.dataPtr(7);
            const Box& obox = ofab.box();
            const int NI=obox.length(0);
            const int NJ=obox.length(1);
            const int NK=obox.length(2);

            const FArrayBox& fab = mf[mfi];
            if ( !(fab.box().contains(obox)))
                BoxLib::Abort();

            const Real* v1x = fab.dataPtr(0);
            const Real* v1y = fab.dataPtr(1);
            const Real* v1z = fab.dataPtr(2);
            const Real* v2x = fab.dataPtr(3);
            const Real* v2y = fab.dataPtr(4);
            const Real* v2z = fab.dataPtr(5);

            for (int k=0; k<NK; k++) {
                for (int j=0; j<NJ; j++) {
                    for (int i=0; i<NI; i++) {

                        const int cell = (k*NJ+j)*NI+i;

                        Real gTmag = v1x[cell]*v1x[cell] + v1y[cell]*v1y[cell] + v1z[cell]*v1z[cell];

                        if (gTmag > epsilon)
                        {
                            Real dot   = v1x[cell]*v2x[cell] + v1y[cell]*v2y[cell] + v1z[cell]*v2z[cell];

                            outx[cell] = v2x[cell] - dot * v1x[cell] / (gTmag + epsilon);
                            outy[cell] = v2y[cell] - dot * v1y[cell] / (gTmag + epsilon);
                            outz[cell] = v2z[cell] - dot * v1z[cell] / (gTmag + epsilon);
                        }
                        else
                        {
                            outx[cell] = 0;
                            outy[cell] = 0;
                            outz[cell] = 0;

                        }

                        outmag[cell] = std::sqrt(outx[cell]*outx[cell] + outy[cell]*outy[cell] + outz[cell]*outz[cell]);

                        if (outmag[cell] <= epsilon)
                        {
                            outx[cell] = 0;
                            outy[cell] = 0;
                            outz[cell] = 0;
                            outmag[cell] = 0;
                        }

                        Real v2mag = std::sqrt(v2x[cell]*v2x[cell] + v2y[cell]*v2y[cell] + v2z[cell]*v2z[cell]);
                        if (v2mag > epsilon)
                        {
                            fx[cell] = outx[cell]/v2mag;
                            fy[cell] = outy[cell]/v2mag;
                            fz[cell] = outz[cell]/v2mag;
                        }
                        else
                        {
                            fx[cell] = 0.;
                            fy[cell] = 0.;
                            fz[cell] = 0.;
                        }
                        fmag[cell] = std::sqrt(fx[cell]*fx[cell] + fy[cell]*fy[cell] + fz[cell]*fz[cell]);


                    } // i
                } // j
            } // k
        } // mfi

        // Fix up fine-fine and periodic
        cross[iLevel].FillBoundary(0,BL_SPACEDIM);
        geoms[iLevel].FillPeriodicBoundary(cross[iLevel],0,BL_SPACEDIM);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& ofab = cross[iLevel][mfi];
            const Real* vecx = ofab.dataPtr(0);
            const Real* vecy = ofab.dataPtr(1);
            const Real* vecz = ofab.dataPtr(2);
            Real* divVec     = ofab.dataPtr(8);

            const Box& box = mfi.validbox();
            const int ni = box.length(0);
            const int nj = box.length(1);
            const int nk = box.length(2);

            for (int k=0; k<nk; k++) {
                for (int j=0; j<nj; j++) {
                    for (int i=0; i<ni; i++) {

                        const int cell = ((k+nGrow)*(nj+2*nGrow)+j+nGrow)*(ni+2*nGrow)+i+nGrow;
                        const int cell_im1 = cell-1;
                        const int cell_ip1 = cell+1;
                        const int cell_jm1 = cell-(ni+2*nGrow);
                        const int cell_jp1 = cell+(ni+2*nGrow);
                        const int cell_km1 = cell-(ni+2*nGrow)*(nj+2*nGrow);
                        const int cell_kp1 = cell+(ni+2*nGrow)*(nj+2*nGrow);

                        Real divx = vecx[cell_ip1] - vecx[cell_im1];
                        Real divy = vecy[cell_jp1] - vecy[cell_jm1];
                        Real divz = vecz[cell_kp1] - vecz[cell_km1];
                        
                        divVec[cell] = -( divx + divy + divz );

                    } // i
                } // j
            } // k
        } // mfi
    } // iLevel

    Array<string> outNames(nCompNEW+comps.size());
    outNames[0] = "parallelFlow_x";
    outNames[1] = "parallelFlow_y";
    outNames[2] = "parallelFlow_z";
    outNames[3] = "parallelFlow_mag";
    outNames[4] = "parallelFlowFrac_x";
    outNames[5] = "parallelFlowFrac_y";
    outNames[6] = "parallelFlowFrac_z";
    outNames[7] = "parallelFlowFrac_mag";
    outNames[8] = "divParallelFlow";
    for (int i=0; i<comps.size(); ++i)
        outNames[nCompNEW+i] = varNames[nCompOLD + i];
    WritePlotFile(cross,outfile,outNames,amrData);
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
