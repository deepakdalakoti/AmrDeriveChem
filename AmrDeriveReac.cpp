#include <sstream>
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
  //  int nReacIn = pp.countval("reacs"); current hack to only get rxns which involve CH2O

    int tot_nReacIn = 268;
    Array<int> tot_Reacs(tot_nReacIn);
    Array<int> Reacs;
    for (int i=0 ; i < tot_nReacIn ; i ++) tot_Reacs[i] = i;
//    pp.getarr("reacs", Reacs);
    int nSpec = cd.numSpecies();
    int sCompX = 0;
    int sCompT = nSpec;
    typedef Array<std::pair <std::string,int> > PairAr;
    Array<PairAr> SpecInR(tot_nReacIn);
    bool hasCH2O;
//    for (int i=0 ; i < nReacIn ; i ++)  {
//       SpecInR[i]=cd.specCoeffsInReactions(Reacs[i]);
//      }
      for (int i=0 ; i < tot_nReacIn ; i ++)  {
       SpecInR[i]=cd.specCoeffsInReactions(tot_Reacs[i]);
           hasCH2O = false;
      
        for (int j=0 ; j < SpecInR[i].size() ; j++)
             if (SpecInR[i][j].first == "CH2O") hasCH2O=true;
          
             if(hasCH2O) {
               Reacs.push_back(i);      
               if(ParallelDescriptor::IOProcessor())
                   for (int k=0 ; k<SpecInR[i].size() ; k++)
                    std::cout << " Reaction  " << i << " " << SpecInR[i][k].first << std:: endl;
                  }

       }
    int nReacIn = Reacs.size();
    if(ParallelDescriptor::IOProcessor()) {
          for  (int i =0 ; i < nReacIn ; i++) std::cout << cd.reactionString(Reacs[i]) << std::endl;
     }
/*  if(ParallelDescriptor::IOProcessor()) {
  for(int j=0 ; j<nReacIn ; j++ ){
     std::cout << "Reaction Number  " << j << std::endl; 
   for(int i=0 ; i < SpecInR[j].size() ; i++ )
   std::cout << SpecInR[j][i].first <<  " " <<SpecInR[j][i].second  << std::endl;
    }
  }
*/
 /*   Array<Real> MW(nSpec) ;
    MW = cd.speciesMolecWt();

    Array<string> varNames(nSpec+1+1);
    for (int i=0; i<nSpec; ++i)
        varNames[sCompX + i] = "Y(" + cd.speciesNames()[i] + ")";
    varNames[sCompT] = "temp";
    varNames[sCompT+1] = "density";

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
        Out.set(iLevel,new MultiFab(amrData.boxArray(iLevel),2*nReacIn,0,Fab_allocate));

        MultiFab mf(amrData.boxArray(iLevel),varNames.size(),0,Fab_allocate);
        amrData.FillVar(mf,iLevel,varNames,destFillComps);
  

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read on level " << iLevel << std::endl;
            FArrayBox H, RR, T_fal,T, Rfw, Rrv,X;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {   const Box box = mfi.validbox();
            FArrayBox& Y = mf[mfi];
            T.resize(box,1);
            T.copy(mf[mfi],nSpec,0,1);
            H.resize(box,nSpec);
            RR.resize(box,nSpec);
            Rfw.resize(box,nReacIn);
            Rrv.resize(box,nReacIn);
            X.resize(box,nSpec);
            FArrayBox& out = Out[iLevel][mfi];
            out.setVal(0.0);
            cd.getHGivenT(H,T,box,0,0); // comes in mass units
            for (int i=0 ; i <nSpec; i++)
            H.mult(MW[i],i,1);
            cd.massFracToMoleFrac(X,Y,box,0,0);
            cd.fwdRevReacRatesGivenXTP(Rfw,Rrv,Reacs,X,T,Patm,box,0,0,0,0);
            out.copy(Rfw,0,0,nReacIn);
            out.minus(Rrv,0,0,nReacIn);
            for  (int i=0 ; i < nReacIn ; i++) {
               for (int j=0 ; j < SpecInR[i].size(); ++j) {
                for (int k=0 ; k < nSpec ; k++ ) {
                   if(SpecInR[i][j].first==cd.speciesNames()[k]) {
                       RR.copy(H,k,0,1);
                       RR.mult(SpecInR[i][j].second,0,1);
                       out.plus(RR,0,nReacIn+i,1);
                           }
                        }

                     }
                   
                       out.mult(out,i,nReacIn+i,1);
                  }
           out.mult(-1.0,nReacIn,nReacIn); 


        }
        if(ParallelDescriptor::IOProcessor())
             std::cout << "Derive Completed on Level " << iLevel << std::endl; 

    }
     std::string ofile("plt_HR");
     pp.query("outfile", ofile);
     Array<std::string> names(2*nReacIn);
     for (int i=0 ; i<nReacIn ; i++)
{
          std::ostringstream oss ;
           oss << "Reaction_Rate_" << i ;
           names[i] = oss.str();
    }
     for (int i=0 ; i<nReacIn ; i++) 
   {    
           std::ostringstream oss;
           oss << "HRR_" << i ;
           names[nReacIn+i] = oss.str();


   }
     AmrData& a = amrData;
     WritePlotFile(Out,ofile,names,a);
*/
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
