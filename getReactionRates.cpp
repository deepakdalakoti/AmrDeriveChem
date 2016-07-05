
#include <iosfwd>
#include <ccse-mpi.H>
#include <SPACE.H>
#include <Array.H>
#include <BLassert.H>
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "VisMF.H"
#include "Utility.H"
#include "ChemDriver.H"
#include "getReactionRates_F.H"
#include "getJPDF_F.H"
#include "preProcessing.H"

typedef ChemDriver::Edge Edge;
typedef std::list<Edge> EdgeList;

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

    ChemDriver cd;
    
    Real Patm=1; pp.query("Patm",Patm);
    std::string infile; pp.get("infile",infile);
    
    
    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
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

    PArray<MultiFab> outState(Nlev,PArrayManage);
    int nComp = 1;
    int nSpec = cd.numSpecies();
    int sCompX = 0;
    int sCompT = nSpec;
    int sCompHeatRelease = nSpec +1;
    Array<string> varNames(nSpec+1+1); // species, temperature, heatrelease
    for (int i=0; i<nSpec; ++i)
       {
          varNames[sCompX + i] = "X(" + cd.speciesNames()[i] + ")";
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
    
    nComp = nSpec + 2; //number of species + temperature + heat release
  
   bool writeout = true;
    
    for (int iLevel=0; iLevel<bas.size(); ++iLevel){
       MultiFab mf(bas[iLevel],varNames.size(),0,Fab_allocate);
       amrData.FillVar(mf,iLevel,varNames,destFillComps);
       
       outState.set(iLevel,new MultiFab(bas[iLevel],nSpec+4*Nreacs+1+1+1,0)); // also store "heat of reaction" in outState. 
       
       // Build volume at this level, use our own dx
       Real vol = 1;
       for (int i=0; i<BL_SPACEDIM; ++i) {
          vol *= amrData.ProbSize()[i] / amrData.ProbDomain()[iLevel].length(i);
       }
      
       
       for (MFIter mfi(mf); mfi.isValid(); ++mfi)
          {
             FArrayBox &fab = mf[mfi];
             const Box& box = mfi.validbox();
             
             FArrayBox Hof(box,nSpec);
             FArrayBox T(box,1);
             T.setVal(298.15);
             FArrayBox totHofProd(box, 1);
             FArrayBox totHofReact(box, 1);
             FArrayBox HoReact(box, 1);
             FArrayBox tmp(box, 1);
             FArrayBox tmpY(box, 1);
             FArrayBox dummy(box, 1);

             Hof.setVal(0.0);
             dummy.setVal(0.0);
             /*get forward, backward and net reaction rates */
             cd.fwdRevReacRatesGivenXTP(outState[iLevel][mfi], outState[iLevel][mfi],reacIds, fab,fab,
                                        Patm,box,sCompX,sCompT,sCompT+1+1,sCompT+1+1+Nreacs);
            
             /* get the heat of formation for each species. */
             cd.getHGivenT(Hof, T, box, 0, 0);
             /* converting J/Kg to J/mol */
             for(int is = 0; is < nSpec; ++is){
                Hof.mult(cd.speciesMolecWt()[is]/1000, is, 1);
             }
             
             cd.moleFracToMassFrac(outState[iLevel][mfi],fab,box,sCompX,sCompX);
           
             
             if (iLevel < bas.size()-1)
                {
                   BoxArray baf = BoxArray(bas[iLevel+1]).coarsen(amrData.RefRatio()[iLevel]);	  
                   std::vector< std::pair<int,Box> > isects = baf.intersections(box);                    
                   for (int ii = 0; ii < isects.size(); ii++)
                      {
                         fab.setVal(0,isects[ii].second,0,nComp);
                         outState[iLevel][mfi].setVal(0,isects[ii].second,sCompT+1+1,Nreacs);
                         outState[iLevel][mfi].setVal(0,isects[ii].second,sCompT+1+1+Nreacs,Nreacs);
                         Hof.setVal(0,isects[ii].second,0,nSpec);
                      }
                }
         //     else if(iLevel==bas.size()-1)
//                 {
//                    IntVect idx(617,895);
//                    if (box.contains(idx)) 
//                       {
//                          std::cout << "SAMPLE" << std::endl;
//                          std::cout << fab(idx,sCompT) << std::endl;
//                          for (int k=0; k<nSpec; ++k){
//                             std::cout << fab(idx,sCompX+k) << std::endl;
//                          }
//                          std::cout << "SAMPLE DONE" << std::endl;                         
//                       }
                      
//                 }
             

   
             FORT_GETSUMRR(box.loVect(),box.hiVect(),
                           outState[iLevel][mfi].dataPtr(sCompT+1+1),    ARLIM(outState[iLevel][mfi].loVect()),  ARLIM(outState[iLevel][mfi].hiVect()),
                           outState[iLevel][mfi].dataPtr(sCompT+1+1+Nreacs),    ARLIM(outState[iLevel][mfi].loVect()),      ARLIM( outState[iLevel][mfi].hiVect()),
                           outState[iLevel][mfi].dataPtr(sCompT+1+1+Nreacs+Nreacs), ARLIM( outState[iLevel][mfi].loVect()), ARLIM( outState[iLevel][mfi].hiVect()), &Nreacs);
             outState[iLevel][mfi].copy(fab, sCompT, sCompT, 1);
             outState[iLevel][mfi].copy(fab, sCompT+1, sCompT+1, 1);
             
             /* calculate the heat of reaction for each reaction. */
             for(int iR = 0; iR < Nreacs; ++iR){
                HoReact.setVal(0.0);
                tmp.setVal(0.0);

                Array<std::pair<std::string,int> > coeffs;
                coeffs = cd.specCoeffsInReactions(iR);
                int size =  coeffs.size();
                
                for(int is = 0; is< size; ++is){
                   for(int js = 0; js < nSpec; ++js){
                      if(coeffs[is].first == cd.speciesNames()[js]){
                         tmp.copy(Hof, js, 0, 1);
                         tmp.mult(coeffs[is].second);
                         HoReact += tmp;
                      }
                   }
                }
                
                HoReact.mult(-1);
                //    /*compute the heat of reaction for this reaction (taking care of sign "+") */
                //                if (writeout) std::cout<<iR<<" "<<*totHofReact.dataPtr(0)<<std::endl;
                
                dummy.copy(outState[iLevel][mfi],sCompT+1+1+2*Nreacs+iR,0,1);
                HoReact *=dummy;
                /* store this heat of reaction in outState*/
                outState[iLevel][mfi].copy(HoReact, 0, sCompT+1+1+3*Nreacs+iR, 1);
                outState[iLevel][mfi].plus(outState[iLevel][mfi],  sCompT+1+1+3*Nreacs+iR,  sCompT+1+1+4*Nreacs, 1);
              
             
             }

             writeout = false;
             
             /*end. outState: Y_i (mass fraction of species), T, HeatRelease, Fwd_rr (nSpec), Bkw_rr(nSpec), rr(nSpec)*/
             /* calculate the heat of reaction for each reaction. */
             
          }

       if (ParallelDescriptor::IOProcessor())
          std::cerr << "Data has been completed on level " << iLevel << std::endl;
     
    }
   
    std::string nameFile ="chem-H-84-reactions-names"; pp.query("nameFile",nameFile);
    std::ifstream istr(nameFile.c_str(),std::ios::in);
    std::list<std::string> rrnames;
    while (istr){
       std::string name;
       istr >> name;
       if (name.size()>0)
          rrnames.push_back(name);
    }
    
    std::string outfile(getFileRoot(infile) + "_reactionRates"); pp.query("outfile",outfile);
    Array<std::string> names(nSpec+3+4*Nreacs);
    for (int i=0; i<nSpec; ++i) {
       names[i] = "Y(" + cd.speciesNames()[i] + ")";
    }
    names[nSpec] = "temp";
    names[nSpec+1] = "HeatRelease";
    int fwdN = nSpec+2;
    std::list<std::string>::const_iterator it=rrnames.begin();
    for (int i= fwdN; i<(fwdN+Nreacs); ++i) {
       names[i] = "F:" + (*it);
       ++it;
    }
    int bwdN = fwdN+Nreacs;
    it=rrnames.begin();
    for (int i= bwdN; i<(bwdN+Nreacs); ++i) {
       names[i] = "B:" + (*it) ;
       ++it;
    }
    int rrN = bwdN+Nreacs;
    it=rrnames.begin();
    for (int i= rrN; i<(rrN+Nreacs); ++i) {
       names[i] = "F-B:" + (*it);
       ++it;
    }
    int hfN = rrN+Nreacs;
    it=rrnames.begin();
    for (int i= hfN; i<(hfN+Nreacs); ++i) {
       names[i] = "HoReact" + (*it);
       ++it;
    }

    names[nSpec+2+4*Nreacs] = "TotHofR";
    
    bool verb=true;
    AmrData& a = amrData;
    writePlotfile(outState, a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),outfile,names,verb);
    
    if (ParallelDescriptor::IOProcessor())
       std::cerr << "Data has been written to a plot file: "<<outfile.c_str()<< std::endl;

   BoxLib::Finalize();
    return 0;
}

