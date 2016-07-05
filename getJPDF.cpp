
#include <Array.H>
#include <FArrayBox.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <preProcessing.H>
#include <DataServices.H>
#include <Utility.H>
#include <getJPDF_F.H>

#ifdef WIN32
static const char* path_sep_str = "\\";
#else
static const char* path_sep_str = "/";
#endif


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
    std::string outfile(getFileRoot(infile) + "_JPDF"); pp.query("outfile",outfile);
   
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);    
    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Box domain = amrData.ProbDomain()[0];
    vector<Real> bbll,bbur;
    if (int nx=pp.countval("bounds"))
    {
        Array<Real> barr;
        pp.getarr("bounds",barr,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(barr.size()==2*d);
        bbll.resize(d);
        bbur.resize(d);
        for (int i=0; i<d; ++i)
        {
            bbll[i] = barr[i];
            bbur[i] = barr[d+i];
        }

        // Find coarse-grid coordinates of bounding box, round outwardly
        for (int i=0; i<BL_SPACEDIM; ++i) {
           const Real dx = amrData.ProbSize()[i] / amrData.ProbDomain()[0].length(i);            
           domain.setSmall(i,std::max(domain.smallEnd()[i], (int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx)/dx)));
           domain.setBig(i,std::min(domain.bigEnd()[i], (int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx)/dx)));
        }
    }
    
    // Build boxarrays for fillvar call
    Box levelDomain = domain;
    Array<BoxArray> bas(Nlev);
    for (int iLevel=0; (iLevel<=finestLevel)&&(bas.size()==Nlev); ++iLevel)
    {
        BoxArray baThisLev = BoxLib::intersect(amrData.boxArray(iLevel),levelDomain);

        if (baThisLev.size() > 0) {
            bas.set(iLevel,baThisLev);
            if (iLevel < finestLevel) {
                levelDomain.refine(amrData.RefRatio()[iLevel]);
            }
        }
        else
        {
            bas.resize(iLevel);
        }
        //std::cout << "lev,ba: " << iLevel << ", " << bas[iLevel] << std::endl;
    }
    
    int nX = 64; pp.get("nX", nX);
    int nY = 64; pp.get("nY", nY);


    Real minX,maxX, minY, maxY; 
    std::string varX = "temp"; pp.query("varX",varX);
    amrData.MinMax(amrData.ProbDomain()[finestLevel], varX, finestLevel, minX, maxX);

    Real deltaX = (maxX - minX)/nX;

    Array<std::string> varYs; pp.queryarr("varYs",varYs);
    Box binBox(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(nX-1, nY-1,0)));
    FArrayBox bins(binBox,1);

    PArray<MultiFab> mfs(Nlev,PArrayManage);
    for (int iLevel=0; iLevel<bas.size(); ++iLevel){
        mfs.set(iLevel, new MultiFab(bas[iLevel],2,0,Fab_allocate));
        MultiFab& mf = mfs[iLevel];
        amrData.FillVar(mf,iLevel,varX,0);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox &fab = mf[mfi];
            const Box& box = mfi.validbox();
            
            if (iLevel < bas.size()-1)
            {
                BoxArray baf = BoxArray(bas[iLevel+1]).coarsen(amrData.RefRatio()[iLevel]);	  
                std::vector< std::pair<int,Box> > isects = baf.intersections(box);                    
                for (int ii = 0; ii < isects.size(); ii++)
                {
                    fab.setVal(0,isects[ii].second,0,1);
                }
            }
        }
    }

    bool first = true;
    for (int i=0; i<varYs.size(); ++i) {

        const std::string& varY = varYs[i];        

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Binning " << varY << " against " << varX << std::endl;
        }

        amrData.MinMax(amrData.ProbDomain()[finestLevel], varY, finestLevel, minY, maxY);
        Real deltaY = Real((maxY-minY)/nY);
    
        bins.setVal(0.0);
        
        for (int iLevel=0; iLevel<bas.size(); ++iLevel){

            MultiFab& mf = mfs[iLevel];
            amrData.FillVar(mf,iLevel,varY,1);
            
            // Build volume at this level, use our own dx
            Real vol = 1;
            for (int i=0; i<BL_SPACEDIM; ++i) {
                vol *= amrData.ProbSize()[i] / amrData.ProbDomain()[iLevel].length(i);
            }
            
            
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                FArrayBox &fab = mf[mfi];
                const Box& box = mfi.validbox();
                
                if (iLevel < bas.size()-1)
                {
                    BoxArray baf = BoxArray(bas[iLevel+1]).coarsen(amrData.RefRatio()[iLevel]);	  
                    std::vector< std::pair<int,Box> > isects = baf.intersections(box);                    
                    for (int ii = 0; ii < isects.size(); ii++)
                    {
                        fab.setVal(0,isects[ii].second,1,1);
                    }
                }                
                
                FORT_GETJPDF(box.loVect(),box.hiVect(),
                             fab.dataPtr(0), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                             fab.dataPtr(1), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                             bins.dataPtr(), ARLIM(bins.loVect()), ARLIM(bins.hiVect()),
                             &minX, &maxX, &deltaX, &nX,
                             &minY, &maxY, &deltaY, &nY, &vol);
                
            }
        }

        int numPts = bins.box().numPts();
        ParallelDescriptor::ReduceRealSum(bins.dataPtr(),numPts);

        Real total = bins.sum(0);
        bins.mult(1/total);
        
        if (ParallelDescriptor::IOProcessor()) {

            // Build out folder if first time through

            if (first) {

                std::cout << "Building out folder: " << outfile << std::endl;

                if (!BoxLib::UtilCreateDirectory(outfile, 0755)) {
                    BoxLib::CreateDirectoryFailed(outfile);
                }
                first=false;
            }


            /*create a dat file to store those min max deltas */
            
            int nDigits = std::log10(varYs.size()+1);
            std::string compFileName = BoxLib::Concatenate(outfile + path_sep_str, i, nDigits);

            std::ofstream osfab(std::string(compFileName + ".fab").c_str());
            bins.writeOn(osfab);
            osfab.close();

            std::ofstream ostxt; ostxt.open(std::string(compFileName + ".dat").c_str());

            ostxt <<"minX = "<<minX<<"\n";
            ostxt <<"maxX = "<<maxX<<"\n";
            ostxt <<"minY = "<<minY<<"\n";
            ostxt <<"maxY = "<<maxY<<"\n";
            ostxt <<"nX = "<<nX<<"\n";
            ostxt <<"nY = "<<nY<<"\n";
            ostxt <<"deltaX = "<<deltaX<<"\n";
            ostxt <<"deltaY = "<<deltaY<<"\n";
            ostxt <<"varX = "<<varX<<"\n";
            ostxt <<"varY = "<<varY<<"\n";
            ostxt <<"ilo = "<<0<<"\n";
            ostxt <<"jlo = "<<0<<"\n";
            ostxt <<"ihi = "<<nX-1<<"\n";
            ostxt <<"jhi = "<<nY-1<<"\n";
            ostxt.close();        

        }
        
    }
    

    BoxLib::Finalize();
    return 0;
}

