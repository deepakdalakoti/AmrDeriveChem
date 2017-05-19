#include "Array.H"
#include "REAL.H"
#include "Box.H"
#include "PArray.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"
#include "AppendToPlotFile.H"

#undef MOLE_FRAC_IN_PLOTFILE
#define MASS_FRAC_IN_PLOTFILE
static
void
WritePlotFile (const PArray<MultiFab>&   mfout,
               const std::string&        ofile,
               const Array<std::string>& names,
               AmrData&                  amrData);


void
appendToPlotFile(AmrData&                  amrData,
                 const PArray<MultiFab>&   mfout,
                 std::string&              oFile,
                 const Array<std::string>& nnames,
                 const std::string&        mfBaseName,
                 bool                      verbose)
{
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    std::string nFileHeader(oFile);
    nFileHeader += "/NewHeader";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    std::ifstream is;

    os.precision(17);
    int finestLevel = mfout.size() - 1;
    int n_var_copy; 
    if (ParallelDescriptor::IOProcessor())
    {
        if (verbose)
        {
            std::cout << "Opening files = " << oFileHeader << " and " << nFileHeader << '\n';
        }
        is.open(oFileHeader.c_str(), std::ios::in|std::ios::binary);
        os.open(nFileHeader.c_str(), std::ios::out|std::ios::binary);
        
        if (os.fail())
            BoxLib::FileOpenFailed(oFileHeader);
        if (is.fail())
            BoxLib::FileOpenFailed(nFileHeader);
        //
        // Start writing plotfile.
        //
        std::string version;
        is >> version;
        os << version << '\n';
        int n_var;
        is >> n_var;
        n_var_copy = n_var;
        os << n_var+nnames.size() << '\n';
        Array<std::string> inames(n_var);
        for (int n = 0; n < n_var; n++) is >> inames[n];
        for (int n = 0; n < n_var; n++) os << inames[n] << '\n';
        for (int n = 0; n < nnames.size(); n++) os << nnames[n] << '\n';
        
        int sdim;
        is >> sdim;
        os << sdim << '\n';
        
        Real time;
        is >> time;
        os << time << '\n';
        
        int oFinestLevel;
        is >> oFinestLevel;
        
        BL_ASSERT(oFinestLevel>=finestLevel);
        os << finestLevel << '\n';
        
        Array<Real> probLo(BL_SPACEDIM);
        for (int i = 0; i < BL_SPACEDIM; i++) is >> probLo[i];
        for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
        os << '\n';
        
        Array<Real> probHi(BL_SPACEDIM);
        for (int i = 0; i < BL_SPACEDIM; i++) is >> probHi[i];
        for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
        os << '\n';
        
        Array<int> refRatio(oFinestLevel);
        for (int i = 0; i < oFinestLevel; i++) is >> refRatio[i];
        for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
        os << '\n';
        
        Array<Box> probDomain(oFinestLevel+1);
        for (int i = 0; i <= oFinestLevel; i++) is >> probDomain[i];
        for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
        os << '\n';
        
        int tmpI;
        for (int i = 0; i <= oFinestLevel; i++) is >> tmpI;
        for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
        os << '\n';
        
        Real dx[BL_SPACEDIM];
        for (int i = 0; i <= oFinestLevel; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                is >> dx[k];
            }
            if (i<=finestLevel)
            {
                for (int k = 0; k < BL_SPACEDIM; k++)
                {
                    os << dx[k] << ' ';
                }
                os << '\n';
            }
        }
        
        int coordSys;
        is >> coordSys;
        os << coordSys << '\n';
        
        int bndry;
        is >> bndry;
        os << bndry << '\n'; // The bndry data width.
    }
    //
    // Write out level by level.
    //
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Write state data.
        //
        int nGrids = amrData.boxArray(iLevel).size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if (ParallelDescriptor::IOProcessor())
        {
            int ilev,ngrd;
            Real time;
            is >> ilev >> ngrd >> time;
            os << ilev << ' ' << ngrd << ' ' << time << '\n';

            int tmpI;
            is >> tmpI;
            os << tmpI << '\n';
            
            Real glocl,gloch;
            for (int i = 0; i < nGrids; ++i)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    is >> glocl >> gloch;
                    os << glocl
                       << ' '
                       << gloch
                       << '\n';
                }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            // NOTE: should already exist!
            //
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
//            if (!BoxLib::UtilCreateDirectory(Level, 0755))
//                BoxLib::CreateDirectoryFailed(Level);
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += '/';
        PathName += mfBaseName;
        
        if (ParallelDescriptor::IOProcessor())
        {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += mfBaseName;
           int idx = 0;
            std::string oldMFname;
           is >> oldMFname;
          std::string c = "L" ;
          std::string test;
          test=oldMFname;
          while (oldMFname[0]==c[0]) {
//            is >> oldMFname;
            os << oldMFname << '\n';
            is >> oldMFname;
//            std::cout << oldMFname << std::endl;         
         
            if(test==oldMFname){
              
//              std::cout << "putback" << oldMFname[0] << std::endl;
             break;    }  
             test = oldMFname;
//            std::cout << test << std::endl;     
            
           }
            is.putback(oldMFname[0]);
            os << RelativePathName << '\n';
        }
        VisMF::Write(mfout[iLevel], PathName, VisMF::OneFilePerCPU);
    }
    
    os.close();
    is.close();
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

  std::string infile; pp.get("infile",infile);
  std::string outfile = infile + std::string("_xi"); pp.query("outfile",outfile);
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(infile, fileType);
  if (!dataServices.AmrDataOk())
    DataServices::Dispatch(DataServices::ExitRequest, NULL);    
  AmrData& amrData = dataServices.AmrDataRef();
  int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel + 1;
  ChemDriver cd;
  int nElts = cd.numElements();
  int nSpec = cd.numSpecies();
  const Array<std::string>& elementNames = cd.elementNames();
  const Array<std::string>& speciesNames = cd.speciesNames();
  Array<Real> speciesMolecWt = cd.speciesMolecWt();
  Array<Real> elementWt = cd.elementAtomicWt();
  int nCompIn = nSpec + 1;
  
 Array<std::string> varNames(nCompIn);
  for (int i=0; i<nSpec; ++i) {
#ifdef MASS_FRAC_IN_PLOTFILE
    varNames[i] = "Y(" + cd.speciesNames()[i] + ")";
#else
    varNames[i] = "X(" + cd.speciesNames()[i] + ")";
#endif
  }

       

  varNames[nSpec] = "temp";
  int iXi = nCompIn;
  int iZ = iXi + 1;
  int iZstar = iZ + nElts;
  int nCompOut = iZstar + 1;

  Array<int> comps(varNames.size());
  for (int i=0; i<comps.size(); ++i) {
    comps[i] = amrData.StateNumber(varNames[i]);
  }

  if (verbose && ParallelDescriptor::IOProcessor())
    for (int i=0; i<varNames.size(); ++i)
      std::cout << "Getting component: " << varNames[i] << std::endl;

  Array<int> destFillComps(varNames.size());
  for (int i=0; i<varNames.size(); ++i)
    destFillComps[i] = i;

  int nC = cd.indexElt("C"); BL_ASSERT(nC >= 0  && nC < nElts);
  int nH = cd.indexElt("H"); BL_ASSERT(nH >= 0  && nH < nElts);
  int nO = cd.indexElt("O"); BL_ASSERT(nO >= 0  && nO < nElts);
  Array<Real> beta(nElts, 0);
  beta[nC] = 2/elementWt[nC];
  beta[nH] = 1/(2*elementWt[nH]);
  beta[nO] = -1/elementWt[nO];

  Array<Array<Real> > mu(nSpec,Array<Real>(nElts));
  for (int spec = 0; spec<nSpec; ++spec) {
    for (int elt = 0; elt<nElts; ++elt) {
      mu[spec][elt] = Real(cd.numberOfElementXinSpeciesY(elementNames[elt],speciesNames[spec]))
       * Real(elementWt[elt]) / Real(speciesMolecWt[spec]);
    }
  }

  Array<Real> X_fu(nSpec,0);
  std::string fuelName = "CH3OCH3"; pp.query("fuelName",fuelName);
  Real Xfu = 0.2; pp.query("Xfu",Xfu);
  BL_ASSERT(Xfu<=1 && Xfu>=0);
    X_fu[cd.index(fuelName)] = Xfu;
  X_fu[cd.index("N2")] = 1 - Xfu;    // currently assumed only fuel and N2 in fuel stream
  Array<Real> Y_fu = cd.moleFracToMassFrac(X_fu);
  Array<Real> Z_fu(nElts,0);

  for (int elt=0; elt<nElts; ++elt) {
    for (int spec=0; spec<nSpec; ++spec) {
      if (mu[spec][elt] != 0) {
        Z_fu[elt] += mu[spec][elt]* Y_fu[spec];

    //    std::cout << "mu_FU["<< speciesNames[spec] << "][" << elementNames[elt]<< "] = " << mu[spec][elt] << std::endl;

      }
    }
  }
  Real Zstar_fu 
    = beta[nC]*Z_fu[cd.indexElt("C")]
    + beta[nH]*Z_fu[cd.indexElt("H")]
    + beta[nO]*Z_fu[cd.indexElt("O")];
  
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Zstar_fu = " << Zstar_fu << std::endl;
  }
  Array<Real> X_ox(nSpec,0);
  X_ox[cd.index("O2")] = 0.15;  // change as per ambient conditions
  X_ox[cd.index("N2")] = 0.85;  
  Array<Real> Y_ox = cd.moleFracToMassFrac(X_ox);
  if(ParallelDescriptor::IOProcessor()) {
  std::cout << "Mass frac of O2 in ox: " << Y_ox[cd.index("O2")] << std::endl;}
  Array<Real> Z_ox(nElts,0);
  for (int elt=0; elt<nElts; ++elt) {
    for (int spec=0; spec<nSpec; ++spec) {
      if (mu[spec][elt] != 0) {
        Z_ox[elt] += mu[spec][elt] * Y_ox[spec];
      }
    }
  }
  Real Zstar_ox 
    = beta[nC]*Z_ox[cd.indexElt("C")]
    + beta[nH]*Z_ox[cd.indexElt("H")]
    + beta[nO]*Z_ox[cd.indexElt("O")];
  
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Zstar_ox = " << Zstar_ox << std::endl;
  }

  // Find stoichiometric value of xi
  int numC = cd.numberOfElementXinSpeciesY("C",fuelName);
  int numH = cd.numberOfElementXinSpeciesY("H",fuelName);
  int numO = cd.numberOfElementXinSpeciesY("O",fuelName);

  Real numOfromCO2 = numC * 2;
  Real numOfromH2O = numH / 2;
  Real numO2stoich = 0.5 * (numOfromCO2 + numOfromH2O - numO);
  Array<Real> X_st(nSpec,0);

  X_st[cd.index(fuelName)] = 1/(1 + numO2stoich*(1 + X_ox[cd.index("N2")]/X_ox[cd.index("O2")]));
  X_st[cd.index("O2")] = numO2stoich * X_st[cd.index(fuelName)];
  X_st[cd.index("N2")] = 1 - X_st[cd.index(fuelName)] - X_st[cd.index("O2")];
  Array<Real> Y_st = cd.moleFracToMassFrac(X_st);
  Array<Real> Z_st(nElts,0);
  for (int elt=0; elt<nElts; ++elt) {
    for (int spec=0; spec<nSpec; ++spec) {
      if (mu[spec][elt] != 0) {
        Z_st[elt] += mu[spec][elt] * Y_st[spec];
      }
    }
  }
  Real Zstar_st 
    = beta[nC]*Z_st[cd.indexElt("C")]
    + beta[nH]*Z_st[cd.indexElt("H")]
    + beta[nO]*Z_st[cd.indexElt("O")];
  
  
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "fuel name " << fuelName << std::endl ;
   std::cout << "Zstar_st = " << Zstar_st << std::endl;
    std::cout << "xi_st = " << (Zstar_st - Zstar_ox)/(Zstar_fu - Zstar_ox) << std::endl;
  }
  PArray<MultiFab> mf(Nlev);
  PArray<MultiFab> Mixfrac(Nlev,PArrayManage);
  std::vector< std::pair<int,Box> > isects;
  FArrayBox tmp, Y;
  for (int iLevel=finestLevel; iLevel>=0; --iLevel)
  { if(ParallelDescriptor::IOProcessor()) std::cout << "Level " << iLevel << std::endl;
    int cr = 1;
    const BoxArray& ba = amrData.boxArray(iLevel);
    mf.set(iLevel,new MultiFab(ba,nCompOut,0,Fab_allocate));
    Mixfrac.set(iLevel,new MultiFab(ba,1,0,Fab_allocate));
    mf[iLevel].setVal(0);
    Mixfrac[iLevel].setVal(0);
    amrData.FillVar(mf[iLevel],iLevel,varNames,destFillComps);

    mf[iLevel].setVal(0,iXi,1);
    for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi)
    {
      const FArrayBox& Spec = mf[iLevel][mfi];
      BoxArray ba(mfi.validbox());

      if (iLevel < finestLevel) {
        cr = amrData.RefRatio()[iLevel];
        const BoxArray baf = BoxArray(amrData.boxArray(iLevel+1)).coarsen(cr);
        ba = BoxLib::complementIn(mfi.validbox(),baf);
      }
      for (int i=0; i<ba.size(); ++i) {
        const Box& box = ba[i];
        FArrayBox& Z = mf[iLevel][mfi];
        Y.resize(box,nSpec);
#ifdef MOLE_FRAC_IN_PLOTFILE
        cd.moleFracToMassFrac(Y,Spec,box,0,0);
#else
        Y.copy(Spec,0,0,nSpec);
#endif  
        tmp.resize(box,1);
        for (int elt=0; elt<nElts; ++elt) {        
          Z.setVal(0,box,iZ+elt,1);
          for (int spec=0; spec<nSpec; ++spec) {
            if (mu[spec][elt] != 0) {
              tmp.copy(Y,spec,0,1);
              tmp.mult(mu[spec][elt]);
              Z.plus(tmp,box,box,0,iZ+elt,1);
            }
          }
          tmp.copy(Z,iZ+elt,0,1);
          tmp.mult(beta[elt],0,1);
          mf[iLevel][mfi].plus(tmp,box,box,0,iZstar,1);
        }
        
        mf[iLevel][mfi].copy(mf[iLevel][mfi],box,iZstar,box,iXi,1);
        mf[iLevel][mfi].plus(-Zstar_ox,box,iXi,1);
        mf[iLevel][mfi].mult(1/(Zstar_fu - Zstar_ox),box,iXi,1);
      }
    }

    if (iLevel < finestLevel) {
      BoxList blC;
      const BoxArray baF = BoxArray(amrData.boxArray(iLevel+1)).coarsen(cr);
      for (int i=0; i<ba.size(); ++i) {
        isects = baF.intersections(ba[i]);
        for (int ii = 0; ii < isects.size(); ii++){
          blC.push_back(isects[ii].second);
        }
      }
      blC.simplify();
      BoxArray baC(blC);
      MultiFab xiC(baC,1,0); xiC.setVal(0);
      MultiFab xiF(BoxArray(baC).refine(cr),1,0);
      xiF.copy(mf[iLevel+1],iXi,0,1); // parallel copy
      for (MFIter mfi(xiC); mfi.isValid(); ++mfi) {
        const Box& vbox = mfi.validbox();
        FArrayBox& cfab = xiC[mfi];
        const FArrayBox& ffab = xiF[mfi];
        for (IntVect civ=vbox.smallEnd(), bend=vbox.bigEnd(); civ<=bend; vbox.next(civ)) {
          Box cbox(civ,civ);
          Box fbox = Box(cbox).refine(cr);
          Real& cval = cfab(civ,0);
          for (IntVect fiv=fbox.smallEnd(), fend=fbox.bigEnd(); fiv<=fend; fbox.next(fiv)) {
            cval += ffab(fiv,0);
          }
          cval *= 1/double(fbox.numPts());
        }
      }
      mf[iLevel].copy(xiC,0,iXi,1); // parallel copy
    }
  }
  for (MFIter mfi(mf[0]); mfi.isValid(); ++mfi) {
    if (mfi.validbox().contains(IntVect(D_DECL(0,0,0)))) {
      IntVect iv(D_DECL(0,0,0));
      Array<Real> Y_test(nSpec);
      for (int n=0; n<nSpec; ++n) {
        Y_test[n] = mf[0][mfi](iv,n);
      }
      Array<Real> Z_test(nElts,0);
      for (int elt=0; elt<nElts; ++elt) {
        for (int spec=0; spec<nSpec; ++spec) {
          Z_test[elt] += Y_test[spec] * mu[spec][elt];
        }
      }
      Real Zstar_test 
          = beta[nC]*Z_test[cd.indexElt("C")]
          + beta[nH]*Z_test[cd.indexElt("H")]
          + beta[nO]*Z_test[cd.indexElt("O")];

      Real xi_test = (Zstar_test - Zstar_ox)/(Zstar_fu - Zstar_ox);
    }
  }
  Array<string> outNames(nCompOut);
  for (int i=0; i<varNames.size(); ++i) {
    outNames[i] = varNames[i];
  }

  outNames[iXi] = "Z";
  for (int i=0; i<nElts; ++i) {
    outNames[iZ+i] = "Z("+elementNames[i]+")";
  }
  outNames[iZstar] = "Zstar";
  for (int iLevel=0; iLevel<=finestLevel; ++iLevel)
 {
  MultiFab::Copy(Mixfrac[iLevel],mf[iLevel], nCompOut-6,0,1,0);
  }
  Array<std::string>  newname(1);
  newname[0] = outNames[nCompOut-6];
  bool append = true ;
  pp.query("append",append);
  if (append) 
  appendToPlotFile(amrData,Mixfrac,infile,newname, "Z",0);
  else 
  WritePlotFile(Mixfrac,outfile,newname,amrData);

  BoxLib::Finalize();
  return 0;
}

       // NOTE: should already exist!            }      is >> sdim;
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
