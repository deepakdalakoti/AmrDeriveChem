
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"
#include "AmrDeriveQPD_F.H"
#include "WritePlotFile.H"

typedef ChemDriver::Edge Edge;
typedef std::list<Edge> EdgeList;
struct ELIcompare
{
    // Order EL iters using the underlying edge
    bool operator()(const EdgeList::const_iterator& lhs, const EdgeList::const_iterator& rhs)
        {
            return *lhs < *rhs;
        }
};

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define INSIDE    1
#define OUTSIDE   0
#define BOUND     2
#define HEATRELEASE     0
#define POLYGON     1

struct TPoint{
   TPoint() : x(.0), y(.0) {};
   TPoint(Real x1, Real y1) : x(x1), y(y1) {};
   
   bool operator==(const TPoint& _right)
      {
         return x == _right.x && y == _right.y;
      };
   
   Real x, y;
};


//horizintal left cross over direction algorithm
//-----------------------------------------------
//  bound   |   value that will be returned only if (p) lies on the bound or vertex
int InsidePolygon(const Array<TPoint>& polygon, int N, TPoint& p, int bound)
{
    //cross points count of x
    int __count = 0;

    //neighbour bound vertices
    TPoint p1, p2;

    //left vertex
    p1 = polygon[0];

    //check all rays
    for(int i = 1; i <= N; ++i)
    {
        //point is an vertex
        if(p == p1) return bound;

        //right vertex
        p2 = polygon[i % N];

        //ray is outside of our interests
        if(p.y < MIN(p1.y, p2.y) || p.y > MAX(p1.y, p2.y))
        {
            //next ray left point
            p1 = p2; continue;
        }

        //ray is crossing over by the algorithm (common part of)
        if(p.y > MIN(p1.y, p2.y) && p.y < MAX(p1.y, p2.y))
        {
            //x is before of ray
            if(p.x <= MAX(p1.x, p2.x))
            {
                //overlies on a horizontal ray
                if(p1.y == p2.y && p.x >= MIN(p1.x, p2.x)) return bound;

                //ray is vertical
                if(p1.x == p2.x)
                {
                    //overlies on a ray
                    if(p1.x == p.x) return bound;
                    //before ray
                    else ++__count;
                }

                //cross point on the left side
                else
                {
                    //cross point of x
                    Real xinters = (p.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;

                    //overlies on a ray
                    if(fabs(p.x - xinters) < __DBL_EPSILON__) return bound;

                    //before ray
                    if(p.x < xinters) ++__count;
                }
            }
        }
        //special case when ray is crossing through the vertex
        else
        {
            //p crossing over p2
            if(p.y == p2.y && p.x <= p2.x)
            {
                //next vertex
                const TPoint& p3 = polygon[(i+1) % N];

                //p.y lies between p1.y & p3.y
                if(p.y >= MIN(p1.y, p3.y) && p.y <= MAX(p1.y, p3.y))
                {
                    ++__count;
                }
                else
                {
                    __count += 2;
                }
            }
        }

        //next ray left point
        p1 = p2;
    }

    //EVEN
    if(__count % 2 == 0) return(OUTSIDE);
    //ODD
    else return(INSIDE);
}

void ZeroFabOutsidePolygon(const Array<TPoint>& polygon,
                           int                  sComp,
                           int                  nComp,
                           const Array<Real>&   dx,
                           FArrayBox&           fab)
{   
    const int Np = polygon.size();
    TPoint p;
    const Box& box = fab.box();
    const IntVect& se = box.smallEnd();
    const IntVect& be = box.bigEnd();
    for (IntVect iv = se; iv <= be; box.next(iv))
    {      
        p.x = (iv[0] + 0.5) * dx[0];
        p.y = (iv[1] + 0.5) * dx[1];
        
        if (InsidePolygon(polygon, Np, p, BOUND) != 1) 
        {
            for(int ir = sComp; ir < nComp; ++ir){
                fab(iv, ir) = 0.0;
            }
        }
    }
}

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

    ChemDriver cd;//("tran.asc.drm19");
//    ChemDriver cd("tran.asc.glarSkel");

    std::string infile; pp.get("infile",infile);

    DataServices::SetBatchMode();
    Amrvis:: FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);    
    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    
    std::string progressName = "HeatRelease"; pp.query("progressName",progressName);
    Real progMin = 0; pp.query("progMin",progMin);
    Real progMax = 100000.; pp.query("progMax",progMax);
    Real lowpercent = 0.0; pp.query("lowpercent", lowpercent);
    Real uppercent = 1.0; pp.query("uppercent", uppercent);
    amrData.MinMax(amrData.ProbDomain()[finestLevel], progressName, finestLevel, progMin, progMax);
    Real lowval = progMax*lowpercent;
    Real upval = progMax*uppercent;
    std::string region = "allburning"; pp.query("region",region);
    Real stdEdgeThick = 20; pp.query("stdEdgeThick",stdEdgeThick);  
    
    Box domain = amrData.ProbDomain()[0];
    std::cout << "domain: " << domain << std::endl;
    vector<int> bbll,bbur;
    if (int nx=pp.countval("bounds"))
    {
        Array<int> barr;
        pp.getarr("bounds",barr,0,nx);
        IntVect sm(D_DECL(barr[0],barr[1],barr[2]));
        IntVect bg(D_DECL(barr[BL_SPACEDIM],barr[BL_SPACEDIM+1],barr[BL_SPACEDIM+2]));
        domain &= Box(sm,bg);
        std::cout << "domain: " << domain << std::endl;
        if (ParallelDescriptor::IOProcessor() && verbose)
          std::cout << "domain: " << domain << std::endl;

          
    }

    // Build boxarrays for fillvar call
    //IntVect sm(D_DECL(0,0,500));
    //IntVect bg(D_DECL(90,90,550));
    //domain &= Box(sm,bg);
    std::cout << "domain: " << domain << std::endl; 
    Box levelDomain = domain;
    amrData.MinMax(domain, progressName, finestLevel, progMin, progMax);
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
        if (ParallelDescriptor::IOProcessor() && verbose)
          std::cout << "lev,ba: " << iLevel << ", " << bas[iLevel] << std::endl;
    }

    Real Patm=1; pp.query("Patm",Patm);

    int nSpec = cd.numSpecies();
    int sCompX = 0;
    int sCompT = nSpec;
    Array<string> varNames(nSpec+1+1); // add heat release too
    for (int i=0; i<nSpec; ++i)
    {
        varNames[sCompX + i] = "Y(" + cd.speciesNames()[i] + ")";
        if (amrData.StateNumber(varNames[sCompX + i]) < 0)
            BoxLib::Abort(std::string("Cannot find " + varNames[sCompX + i]).c_str());
    }
    varNames[sCompT] = "temp";
    varNames[sCompT+1] = "HeatRelease";

    Array<int> destFillComps(varNames.size());
    for (int i=0; i<varNames.size(); ++i)
        destFillComps[i] = i;
    
    int Nreacs = cd.numReactions();
    Array<int> reacIds(Nreacs);
    for (int i=0; i<Nreacs; ++i)
    {
        reacIds[i] = i;
    }
    Array<Real> Fsum(Nreacs,0);
    Array<Real> Rsum(Nreacs,0);
    FArrayBox Fwd, Rev, maskFab;

    if(BL_SPACEDIM == 2) {
    int nv = pp.countval("poly");
    if (nv % 2 != 0)
       BoxLib::Abort("Need even number of values for poly");
    const int Np = nv/2;
    Array<Real> poly(nv); pp.getarr("poly",poly);
    Array<TPoint> polygon(Np); 
    for (int ip = 0; ip < Np; ++ip){
       polygon[ip] = TPoint(poly[2*ip],poly[2*ip+1]);
    }
    }

    std::string outfile_check(""); pp.query("outfile_check",outfile_check);
    PArray<MultiFab> state(Nlev,PArrayManage);

    Array<Real> DX(BL_SPACEDIM); 
    long totalFineCellsInPoly = 0;
    for (int iLevel=0; iLevel<bas.size(); ++iLevel){
        MultiFab mf(bas[iLevel],varNames.size(),0,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);
        
        if (outfile_check!="") {
            state.set(iLevel,new MultiFab(bas[iLevel],varNames.size(),0));
            state[iLevel].copy(mf);
        }
        
        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;
        
        // Build volume at this level, use our own dx
        Real vol = 1;
        for (int i=0; i<BL_SPACEDIM; ++i) {
            vol *= amrData.ProbSize()[i] / amrData.ProbDomain()[iLevel].length(i);
            DX[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[iLevel].length(i);
        }

        long totalCellsThisLevel = 0;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const FArrayBox& Y = mf[mfi];
            const FArrayBox& T = mf[mfi];
            const FArrayBox& HR = mf[mfi];
            FArrayBox& X = mf[mfi];
            
            const Box& box = mfi.validbox();
            Fwd.resize(box,Nreacs);
            Rev.resize(box,Nreacs);
            maskFab.resize(box,1);
            maskFab.setVal(1.);
            
            cd.massFracToMoleFrac(X,Y,box,0,0);
            cd.fwdRevReacRatesGivenXTP(Fwd,Rev,reacIds,X,T,Patm,box,sCompX,sCompT,0,0);
            
            // Zero out covered data
            if (iLevel < bas.size()-1){
                BoxArray baf = BoxArray(bas[iLevel+1]).coarsen(amrData.RefRatio()[iLevel]);	  
                std::vector< std::pair<int,Box> > isects = baf.intersections(box);                    
                for (int ii = 0; ii < isects.size(); ii++){
                    maskFab.setVal(0,isects[ii].second,0,1);
                    if (outfile_check!="")
                        state[iLevel][mfi].setVal(0,isects[ii].second,0,state[iLevel].nComp());
                }
            }
            
#if     HEATRELEASE 
            /* if heat release not in the specified region, then zero the fwd and rev reactions. */ 
            FORT_CONDITIONEDRR(box.loVect(),box.hiVect(),
                               HR.dataPtr(sCompT+1),ARLIM(HR.loVect()),ARLIM(HR.hiVect()),
                               Fwd.dataPtr(),ARLIM(Fwd.loVect()),ARLIM(Fwd.hiVect()),
                               Rev.dataPtr(),ARLIM(Rev.loVect()),ARLIM(Rev.hiVect()),
                               &lowval, &upval, &Nreacs);
#endif         
/*            
#if      POLYGON       
            ZeroFabOutsidePolygon(polygon,0,1,DX,maskFab);
            if (outfile_check!="")
                ZeroFabOutsidePolygon(polygon,0,state[iLevel].nComp(),DX,state[iLevel][mfi]);
#endif*/

            totalCellsThisLevel += maskFab.sum(0);
            for (int i=0; i<Nreacs; ++i)
            {
                Fwd.mult(maskFab,0,i,1);
                Rev.mult(maskFab,0,i,1);
                
                Fsum[i] += Fwd.sum(i);
                Rsum[i] += Rev.sum(i);
            }
        }
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Number of cells in polygon at level " << iLevel << " " << totalCellsThisLevel << std::endl;

        totalFineCellsInPoly += totalCellsThisLevel;

        if (iLevel<finestLevel) {
            int ratio = amrData.RefRatio()[iLevel];
            totalFineCellsInPoly *= ratio*ratio;
        }
    }
    ParallelDescriptor::ReduceLongSum(totalFineCellsInPoly,ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Number of fine cells in polygon " << totalFineCellsInPoly << std::endl;
    
    if (outfile_check!="")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Writing new data to " << outfile_check << std::endl;

        const AmrData& a = amrData;
        WritePlotfile("NavierStokes-V1.1",state,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                      a.DxLevel(),a.CoordSys(),outfile_check,varNames,0);
    }

    ParallelDescriptor::ReduceRealSum(Fsum.dataPtr(),Nreacs,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(Rsum.dataPtr(),Nreacs,ParallelDescriptor::IOProcessorNumber());

    std::string QPDatom="C"; pp.query("QPDatom",QPDatom);
    std::string QPDlabel = infile; pp.query("QPDlabel",QPDlabel);
    std::string QPDfileName = infile + "_" + region +"_QPD.dat"; pp.query("QPDfileName",QPDfileName);

    if (ParallelDescriptor::IOProcessor())
    {
        std::ofstream osfr(QPDfileName.c_str());
        osfr << QPDlabel << std::endl;
        const Array<std::string>& names = cd.speciesNames();
        for (int i=0; i<names.size(); ++i)
            osfr << names[i] << " ";
        osfr << std::endl;
        
        EdgeList edges = cd.getEdges(QPDatom);
        std::cout<<"\n total edges "<<edges.size()<<std::endl;
        std::map<EdgeList::const_iterator,Real,ELIcompare> F, R;
        Real normVal = 0;
        Real normVal_temp = 0;
        
        for (EdgeList::const_iterator it = edges.begin(); it!=edges.end(); ++it)
        {
            const Array<std::pair<int,Real> > RWL=it->rwl();
            for (int i=0; i<RWL.size(); ++i) {

                F[it] += Fsum[ RWL[i].first ]*RWL[i].second;
                R[it] += Rsum[ RWL[i].first ]*RWL[i].second;

            }
            
            if (it->touchesSp("NC12H26"))// and it->touchesSp("CH3"))
            {
               normVal_temp = stdEdgeThick/(F[it]-R[it]); // Normalize to CH4 destruction on CH4->CH3 edge
               if (it->right()=="NC12H26")
                   normVal_temp *= -1;
               normVal += normVal_temp;
            }
        }
        std::cout << "NormVal: " << normVal << std::endl;

        for (EdgeList::const_iterator it = edges.begin(); it!=edges.end(); ++it)
        {
            if (normVal!=0)
            {
                F[it] *= normVal;
                R[it] *= normVal;
            }
            osfr << it->left() << " " << it->right() << " " << F[it] << " " << -R[it] << " " << '\n'; 
        }
        for (EdgeList::const_iterator it = edges.begin(); it!=edges.end(); ++it)
        {
            if (it->touchesSp("NC12H26"))
            {
                std::cout << *it << std::endl;
                std::map<std::string,Real> edgeContrib;
                int rSgn = ( it->left()=="NC12H26"  ?  -1  :  +1);
                const Array<std::pair<int,Real> > RWL=it->rwl();
                for (int i=0; i<RWL.size(); ++i) {
                    int rxn = RWL[i].first;
                    Real wgt = RWL[i].second;
                    Array<std::pair<std::string,int> > specCoefs = cd.specCoeffsInReactions(rxn);

                    // Find name of reaction partner(s), NP = no partner
                    int thisSgn;
                    for (int j=0; j<specCoefs.size(); ++j)
                    {
                        if (specCoefs[j].first == "NC12H26")
                            thisSgn = specCoefs[j].second;
                    }
                    std::string partnerName = "";
                    for (int j=0; j<specCoefs.size(); ++j)
                    {
                        const std::string& sp = specCoefs[j].first;
                        if (sp != "NC12H26"  &&  thisSgn*specCoefs[j].second > 0)
                            partnerName = ( partnerName != ""  ?  partnerName + "+" + sp : sp);
                    }
                    if (partnerName=="")
                        partnerName="NP";

                    edgeContrib[partnerName] += wgt*(Fsum[rxn]-Rsum[rxn])*normVal/stdEdgeThick;
                }

                Real sump=0, sumn=0;
                for (std::map<std::string,Real>::const_iterator cit=edgeContrib.begin(); cit!=edgeContrib.end(); ++cit)
                {
                    std::cout << "   partner: " << cit->first << " " << cit->second << std::endl;
                    if (cit->second > 0.)
                        sump += cit->second;
                    else
                        sumn += cit->second;
                }
                std::cout << "     sum +ve,-ve: " << sump << " " << sumn << std::endl;
            }
        }
    }

    BoxLib::Finalize();
    return 0;
}

