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

#ifndef WIN32
#include <unistd.h>
#endif
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "FArrayBox.H"
#include "Geometry.H"
#include "preProcessing.H"

#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include "compEnergyTerms_F.H"

using std::cout;
using std::cerr;
using std::endl;



#define TOTALENTHALPY 0
#define SENSIBLEENTHALPY 1


int main (int   argc, char* argv[]){
   BoxLib::Initialize(argc,argv);
   
   if (argc < 2)
      print_usage(argc,argv);
   ParmParse pp;
     
   if (pp.contains("help"))
      print_usage(argc,argv);
     
   int verbose=0; pp.query("verbose",verbose);
   if (verbose>1) AmrData::SetVerbose(true);
     
   std::string plotFileName; pp.get("infile",plotFileName);
   std::string outfile(getFileRoot(plotFileName) + "_energyTerms"); pp.query("outfile",outfile);
   
   //std::string TransportFile="tran.asc.chem-H"; pp.query("TransportFile",TransportFile);
   std::string TransportFile="tran.asc.drm19"; pp.query("TransportFile",TransportFile);
   
   ChemDriver cd(TransportFile);
  
   DataServices::SetBatchMode();
   FileType fileType(NEWPLT);
  
   DataServices dataServices(plotFileName, fileType);
   if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
   }
   AmrData& amrData = dataServices.AmrDataRef();
   int idXin = -1;
   int idTin = -1;
   int idRhoin = -1;
   int idVxin = -1;
   int idVyin = -1;
   int icount = 0;
   

   const Array<std::string>& plotVarNames = amrData.PlotVarNames();
   const std::string spName = "X(" + cd.speciesNames()[0] + ")";
   const int nSpec = cd.numSpecies();
   for (int i=0; i<plotVarNames.size(); ++i)
      {
         if (plotVarNames[i] == spName) idXin = i;
         if (plotVarNames[i] == "temp") idTin = i;
         if (plotVarNames[i] == "density"){
            idRhoin = i;
            icount++;
         }
         if (plotVarNames[i] == "x_velocity"){
            idVxin = i;
            icount++;
         }
         if (plotVarNames[i] == "y_velocity"){
            idVyin = i;
            icount++;
         }
         
      }
   if (ParallelDescriptor::IOProcessor() && (idXin<0 || idTin<0) )
      cerr << "Cannot find required data in pltfile" << endl;
   
   
   int finestLevel = amrData.FinestLevel();
   pp.query("finestLevel",finestLevel);
   int Nlev = finestLevel + 1;
   const int idXst = 0;
   const int idTst = nSpec;
   const int nCompIn = idTst + 1 + icount;
   
   Array<int> destFillComps(nCompIn);
   for (int i=0; i<nCompIn; ++i)
      destFillComps[i] = i;
   Array<std::string> inVarNames(nCompIn);
   for (int i=0; i<nSpec; ++i){
      inVarNames[idXst+i] = plotVarNames[idXin+i];
   }
   inVarNames[idTst] = plotVarNames[idTin];
   inVarNames[idTst+1] = plotVarNames[idRhoin];
   inVarNames[idTst+2] = plotVarNames[idVxin];
   inVarNames[idTst+3] = plotVarNames[idVyin];
      
   PArray<MultiFab> outState(Nlev,PArrayManage);
  
   PArray<Geometry> geoms(Nlev,PArrayManage);
   const int nGrow = 1;
   
   for (int lev=0; lev<Nlev; ++lev){
      const BoxArray ba = amrData.boxArray(lev);
      MultiFab inState(ba,nCompIn,nGrow);
      outState.set(lev,new MultiFab(ba,2*(nSpec+1)+2,0));


      
      Array<Real> dx(BL_SPACEDIM);
      
      for (int i=0; i<BL_SPACEDIM; ++i){
         dx[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);
      }
      
      if (ParallelDescriptor::IOProcessor())
         cerr << "Reading data for level " << lev << endl;
      
      amrData.FillVar(inState,lev,inVarNames,destFillComps);
      for (int i=0; i<inVarNames.size(); ++i){
         amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));
      }
      
      if (ParallelDescriptor::IOProcessor())
         cerr << "Data has been read for level " << lev << endl;
      
      const bool do_corners = false;
      geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));
      
      // Extrap grow cells, as a first guess
      const Box& dbox = amrData.ProbDomain()[lev];
      for (MFIter mfi(inState); mfi.isValid(); ++mfi){
         FArrayBox& s = inState[mfi];
         const Box& box = mfi.validbox();
         FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                       s.dataPtr(),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                       &nCompIn);
      }
          
    
      // Fix up fine-fine and periodic
      inState.FillBoundary(0,nCompIn);
      geoms[lev].FillPeriodicBoundary(inState,0,nCompIn,do_corners);
   
      FArrayBox Y,H,RhoD,Ht,Cpmix;
      for (MFIter mfi(inState); mfi.isValid(); ++mfi)
         {
            const Box& vbox = mfi.validbox();
            const Box gbox = Box(vbox).grow(nGrow);
            H.resize(gbox,nSpec);
            Y.resize(gbox,nSpec);
            Ht.resize(gbox,1);
            mixMolecW.resize(gbox,1);
            Cpmix.resize(gbox,1);
          
            RhoD.resize(gbox,nSpec+1); // last comp will hold lambda
            const FArrayBox& X = inState[mfi];
            const FArrayBox& T = inState[mfi];
            const FArrayBox& Rho = inState[mfi];
            const FArrayBox& Vx = inState[mfi];
            const FArrayBox& Vy = inState[mfi];
            
            FArrayBox& DTerms = outState[lev][mfi];
            Ht.setVal(0.);
            DTerms.setVal(0.);   
            Cpmix.setVal(0.); 
            RhoD.setVal(0.);   
            mixMolecW.setVal(0.);
         
            // Compute some stuff on grown box
            cd.moleFracToMassFrac(Y,X,gbox,idXst,0);
            cd.getHGivenT(H,T,gbox,idTst,0);

            cd.getHmixGivenTY(Ht,T,Y,gbox,idTst,0,0);
            cd.getCpmixGivenTY(Cpmix,T,Y,gbox,idTst,0,0);
         
            /* converting J/Kg to J/mol */
            for(int is = 0; is < nSpec; ++is){
               H.mult(cd.speciesMolecWt()[is]/1000, is, 1);
               
            }
         
            const Real Patm = 1.0; 
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

#if  TOTALENTHALPY 
            FORT_COMPCONVECTERMS(vbox.loVect(),vbox.hiVect(),
                                 Rho.dataPtr(idTst+1),    ARLIM(Rho.loVect()),      ARLIM(Rho.hiVect()),
                                 Y.dataPtr(),         ARLIM(Y.loVect()),  ARLIM(Y.hiVect()),
                                 H.dataPtr(),         ARLIM(H.loVect()),  ARLIM(H.hiVect()),
                                 Ht.dataPtr(),         ARLIM(Ht.loVect()),  ARLIM(Ht.hiVect()),
                                 Vx.dataPtr(idTst+2),        ARLIM(Vx.loVect()),     ARLIM(Vx.hiVect()),
                                 Vy.dataPtr(idTst+3),    ARLIM(Vy.loVect()),      ARLIM(Vy.hiVect()),
                                 DTerms.dataPtr(nSpec+3),    ARLIM(DTerms.loVect()), ARLIM(DTerms.hiVect()),
                                 dx.dataPtr(), &nSpec);
            
            FORT_COMPDIFFTERMS(vbox.loVect(),vbox.hiVect(),
                               T.dataPtr(idTst),    ARLIM(T.loVect()),      ARLIM(T.hiVect()),
                               Y.dataPtr(),         ARLIM(Y.loVect()),      ARLIM(Y.hiVect()),
                               H.dataPtr(),         ARLIM(H.loVect()),      ARLIM(H.hiVect()),
                               RhoD.dataPtr(),      ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               RhoD.dataPtr(nSpec), ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               DTerms.dataPtr(0),    ARLIM(DTerms.loVect()) ,ARLIM(DTerms.hiVect()),
                               dx.dataPtr(), &nSpec);
#endif
            
#if SENSIBLEENTHALPY
     

            FORT_COMPCDDRTERMS(vbox.loVect(),vbox.hiVect(),
                               Rho.dataPtr(idTst+1),    ARLIM(Rho.loVect()),      ARLIM(Rho.hiVect()),
                               Y.dataPtr(),         ARLIM(Y.loVect()),      ARLIM(Y.hiVect()),
                               T.dataPtr(idTst),    ARLIM(T.loVect()),      ARLIM(T.hiVect()),
                               Cpmix.dataPtr(),    ARLIM(Cpmix.loVect()),      ARLIM(Cpmix.hiVect()),
                               H.dataPtr(),         ARLIM(H.loVect()),      ARLIM(H.hiVect()),
                               Ht.dataPtr(),         ARLIM(Ht.loVect()),  ARLIM(Ht.hiVect()),
                               Vx.dataPtr(idTst+2),        ARLIM(Vx.loVect()),     ARLIM(Vx.hiVect()),
                               Vy.dataPtr(idTst+3),    ARLIM(Vy.loVect()),      ARLIM(Vy.hiVect()),
                               RhoD.dataPtr(),      ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               RhoD.dataPtr(nSpec), ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               DTerms.dataPtr(),    ARLIM(DTerms.loVect()) ,ARLIM(DTerms.hiVect()),
                               dx.dataPtr(), &nSpec);
         
#endif
       
            
         }        
      
      if (ParallelDescriptor::IOProcessor())
         cerr << "Derive completed on level " << lev << endl;
      
   }

   if (ParallelDescriptor::IOProcessor())
      cout << "Writing diffusion terms data to " << outfile << endl;
   
   Array<std::string> names(2*(nSpec+1)+2);
   names[2*nSpec] =  "Div(rho.Cp.T.V)";
   names[2*nSpec+1] = "Div(lambda.Grad(T))";
   names[2*nSpec+2] = "sumDiv(D_i.rho.Cp.T.Grad(Y_i))";
   names[2*nSpec+3] = "rhoCp_mix";
   
   for (int i=0; i<nSpec; ++i) {
      names[i] = "Div(rho.D_" + cd.speciesNames()[i] + ".Grad(Y_" + cd.speciesNames()[i] + "))";
      names[nSpec+i] = "Div(rho.V.Y_" + cd.speciesNames()[i] + "))";
     
   }
   bool verb=true;
   AmrData& a = amrData;
   writePlotfile(outState, a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                 a.DxLevel(),a.CoordSys(),outfile,names,verb);
   
  
  

   BoxLib::Finalize();
   return 0;
}
