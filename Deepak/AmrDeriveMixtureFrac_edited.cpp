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
 // std::cout << "checkpoint1" << std::endl;
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(infile, fileType);
  if (!dataServices.AmrDataOk())
    DataServices::Dispatch(DataServices::ExitRequest, NULL);    
  AmrData& amrData = dataServices.AmrDataRef();
  int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel + 1;
// std::cout <<"Checkpoint2" << std::endl;
//  std::cout << "above Chemdriver" << std::endl;
  ChemDriver cd;
//  std::cout << "below chemdriver" << std::endl;
  int nElts = cd.numElements();
//  std::cout << "below reference first to " << std::endl;
  int nSpec = cd.numSpecies();
  const Array<std::string>& elementNames = cd.elementNames();
  const Array<std::string>& speciesNames = cd.speciesNames();
//  std::cout << "above Molecwt"<< std::endl;
 // Array<Real> wt =cd.elementAtomicWt();
  Array<Real> speciesMolecWt = cd.speciesMolecWt();
  Array<Real> elementWt = cd.elementAtomicWt();
//  std::cout << "Deepak " << std::endl;
//  for(int i=0; i <nElts ; ++i) {
//     elementWt[i] = i*2+1;
//}
//  std::cout << "hehe" << cd.elementAtomicWt(0) <<std::endl;
  int nCompIn = nSpec + 1;
//  std::cout << "above string " << nCompIn<<std::endl;
  
 Array<std::string> varNames(nCompIn);
// std::cout << "just above llop" << std::endl; 
  for (int i=0; i<nSpec; ++i) {
#ifdef MASS_FRAC_IN_PLOTFILE
    varNames[i] = "Y(" + cd.speciesNames()[i] + ")";
//    std::cout << varNames[i] << std::endl;
#else
    varNames[i] = "X(" + cd.speciesNames()[i] + ")";
#endif
  }

  //    elementWt = cd.elementAtomicWt();
       

//  std::cout << "above temp" << std::endl; 
  varNames[nSpec] = "temp";
//  std::cout << "below varNames " << std::endl;
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
//  std::cout << nC << nH << nO <<std::endl;
  //Real Carbon = cd.elementAtomicWt(0);
  Array<Real> beta(nElts, 0);
  beta[nC] = 2/elementWt[nC];
//  std::cout << "atomic wt C " << elementWt[nC] << std::endl;
//  std::cout << "atomic wt H " << elementWt[nH];
//  std::cout << "atomic wt O " << elementWt[nO];
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
  X_fu[cd.index("N2")] = 1 - Xfu;
  Array<Real> Y_fu = cd.moleFracToMassFrac(X_fu);
//  std::cout << "Mass frac of fuel: " << Y_fu[cd.index(fuelName)] << std::endl;
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
 //  std::cout << "checkpoint 3 " << std::endl;
  Array<Real> X_ox(nSpec,0);
  X_ox[cd.index("O2")] = 0.21;
  X_ox[cd.index("N2")] = 0.79;
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

  X_st[cd.index(fuelName)] = 1/(1 + numO2stoich*(1 + 0.79/0.21));
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
  PArray<MultiFab> Mixfrac(Nlev);
  std::vector< std::pair<int,Box> > isects;
  FArrayBox tmp, Y;
  for (int iLevel=finestLevel; iLevel>=0; --iLevel)
  {
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
  //    std::cout << "Stuff at ll corner: " << std::endl;
      IntVect iv(D_DECL(0,0,0));
  //    std::cout << "MF: " << mf[0][mfi](iv,iXi) << std::endl;
 //     std::cout << "Y: " << std::endl;
      Array<Real> Y_test(nSpec);
      for (int n=0; n<nSpec; ++n) {
//        std::cout << "  " << speciesNames[n] << " " << mf[0][mfi](iv,n) << " " << Y_ox[n] << std::endl;
        Y_test[n] = mf[0][mfi](iv,n);
      }
      Array<Real> Z_test(nElts,0);
      for (int elt=0; elt<nElts; ++elt) {
        for (int spec=0; spec<nSpec; ++spec) {
          Z_test[elt] += Y_test[spec] * mu[spec][elt];
        }
//        std::cout << "Z_test[" << elementNames[elt] << "] = " << Z_test[elt] << std::endl;
      }
      Real Zstar_test 
          = beta[nC]*Z_test[cd.indexElt("C")]
          + beta[nH]*Z_test[cd.indexElt("H")]
          + beta[nO]*Z_test[cd.indexElt("O")];
//      std::cout << "Zstar_test = " << Zstar_test << std::endl;

      Real xi_test = (Zstar_test - Zstar_ox)/(Zstar_fu - Zstar_ox);
//      std::cout << "xi_test = " << xi_test << std::endl;
    }
  }

  Array<string> outNames(nCompOut);
  for (int i=0; i<varNames.size(); ++i) {
    outNames[i] = varNames[i];
  }

  outNames[iXi] = "Z";
  for (int i=0; i<nElts; ++i) {
    outNames[iZ+i] = "Z("+elementNames[i]+")";
 //   std::cout << elementWt[i] <<std::endl;
  }
  outNames[iZstar] = "Zstar";
  for (int iLevel=0; iLevel<=finestLevel; ++iLevel)
 {
  MultiFab::Copy(Mixfrac[iLevel],mf[iLevel], nCompOut-6,0,1,0);
  }
  Array<std::string>  newname(1);
  newname[0] = outNames[nCompOut-6];
//  WritePlotFile(Mixfrac,outfile,newname,amrData); 
//  std::string newout = "./";
  AppendToPlotFile(amrData,Mixfrac,infile,newname, "Z","NewHeader",1);

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
