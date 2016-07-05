#include <iostream>
#include <cstdio>
using std::cout;
using std::endl;
using std::string;

#include "Utility.H"
#include "ParallelDescriptor.H"
#include "ChemDriver.H"
#include "ParmParse.H"

#ifdef USE_ARRAYVIEW
#include "ArrayView.H"
#endif

const bool    verbose_DEF = true;

void print_state(Real                 T,
		 Real                 Patm,
		 const Array<Real>&   massFrac,
		 const ChemDriver& cd)
{
    int Nspecies = massFrac.size();
    const Array<std::string>& names = cd.speciesNames();
    Array<Real> moleFrac = cd.massFracToMoleFrac(massFrac);
    
    printf("\t             |     Mass      |      Mole\n");
    printf("\t   Species   |   Fraction    |    Fraction\n");
    printf("\t             |               |            \n");
    printf("\t-------------|---------------|--------------\n");
    Real massFracTot, moleFracTot;
    massFracTot = moleFracTot = 0;
    IntVect iv = IntVect::TheZeroVector();
    FArrayBox junk(Box(iv,iv),Nspecies+2);
    for (int i=0; i<Nspecies; i++)
    {
        if (1 || massFrac[i] > 0.0)
        {
            printf("\t%9s    |  %10.4e   |   %10.4e\n",
                   names[i].c_str(),
                   massFrac[i],
                   moleFrac[i]);
	    massFracTot += massFrac[i];
	    moleFracTot += moleFrac[i];
        }
        junk(iv,i) = massFrac[i];
    }
    junk(iv,Nspecies) = T;
    cd.getHmixGivenTY(junk,junk,junk,Box(iv,iv),Nspecies,0,Nspecies+1);
    printf("\t-------------|---------------|--------------\n");
    printf("\t     total   |  %10.4e   |   %10.4e\n", massFracTot, moleFracTot);
    printf("\n\tTemperature(K): %6.2f",T);
    printf("\n\tPressure (atm): %6.2f\n",Patm);
    printf("\n\tHmix (J/kg): %6.2f\n\n",junk(iv,Nspecies+1));
}

void print_rstate(Real                 Temp,
                  Real                 RhoH,
	 	  Real                 Patm,
		  const Array<Real>&   rhoY,
		  const ChemDriver&    cd)
{
    int Nspecies = rhoY.size();
    const Array<std::string>& names = cd.speciesNames();

    Real rho = 0;
    for (int i=0; i<Nspecies; ++i) {
      rho += rhoY[i];
    }
    Array<Real> massFrac(Nspecies);
    for (int i=0; i<Nspecies; ++i) {
      massFrac[i] = rhoY[i] / rho;
    }
    Array<Real> moleFrac = cd.massFracToMoleFrac(massFrac);
    
    printf("\t             |     Mass      |      Mole\n");
    printf("\t   Species   |   Fraction    |    Fraction\n");
    printf("\t             |               |            \n");
    printf("\t-------------|---------------|--------------\n");
    Real massFracTot, moleFracTot;
    massFracTot = moleFracTot = 0;
    IntVect iv = IntVect::TheZeroVector();
    FArrayBox junk(Box(iv,iv),Nspecies+2);
    for (int i=0; i<Nspecies; i++)
    {
        if (1 || massFrac[i] > 0.0)
        {
            printf("\t%9s    |  %10.4e   |   %10.4e\n",
                   names[i].c_str(),
                   massFrac[i],
                   moleFrac[i]);
	    massFracTot += massFrac[i];
	    moleFracTot += moleFrac[i];
        }
        junk(iv,i) = massFrac[i];
    }
    junk(iv,Nspecies) = RhoH / rho;
    junk(iv,Nspecies+1) = Temp;

    //Real errMAX = 1.e-6;
    //cd.getTGivenHY(junk,junk,junk,junk.box(),Nspecies,0,Nspecies+1,errMAX);
    printf("\t-------------|---------------|--------------\n");
    printf("\t     total   |  %10.4e   |   %10.4e\n", massFracTot, moleFracTot);
    printf("\n\tTemperature(K): %6.2f",Temp);
    printf("\n\tPressure (atm): %6.2f\n",Patm);
    printf("\n\tHmix (J/kg): %6.2f\n\n",junk(iv,Nspecies));
}


int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    // Parse command line
    ParmParse pp;

    int vin=(verbose_DEF?1:0); pp.query("verbose",vin);
    bool verbose = (vin==1 ? true : false);
    
    ChemDriver cd;
    if (verbose)
        cd.set_verbose_vode();

    const Array<std::string>& names = cd.speciesNames();
    Array<Real> mwt = cd.speciesMolecWt();

    const int nSpec = cd.numSpecies();

    // Set all locations to same IC
    const int nComp = nSpec + 4;
    const int sCompY = 2;
    const int sCompT = nSpec + sCompY;
    const IntVect probe = IntVect::TheZeroVector();
    int npts = 1; pp.query("npts",npts);
    const Box box(-(npts-1)*BoxLib::BASISV(0),probe);
    FArrayBox ostate(box,nComp);
    FArrayBox nstate(box,nComp);
    ostate.setVal(0.);
    nstate.setVal(0.);

    int sCompR = 0;
    int sCompH = 1;

#if 0
    ostate.setVal(.4,box,sCompY+cd.index("CH4"),1);
    ostate.setVal(.4,box,sCompY+cd.index("O2"),1);
#else
    // Read in from file
    const Array<std::string>& speciesNames = cd.speciesNames();
    ostate.setVal(0,ostate.box(),sCompY,nSpec);

    for (int i=0; i<speciesNames.size(); ++i)
    {
        Real val = 0;
        if (pp.countval(speciesNames[i].c_str())) {
          pp.query(speciesNames[i].c_str(),val);
          ostate.setVal(val,box,sCompY+i,1);
        }
    }
#endif

    bool do_sdc = false; pp.query("do_sdc",do_sdc);
    if (!do_sdc) {
      Real sum=0.0;
      for (int i=0; i<nSpec; i++)
        if (i!=cd.index("N2"))
          sum += ostate(probe,sCompY+i);
      BL_ASSERT(sum<=1.0 && sum>=0.0);
      ostate.setVal(1.0-sum,box,sCompY+cd.index("N2"),1);
    }
    
    // Set Temp, pressure and integration time
    Real Temp = 1500; pp.query("T",Temp);
    Real Patm = 1; pp.query("Patm",Patm);
    Real dt = 3.e-5; pp.query("dt",dt);
    ostate.setVal(Temp,sCompT);

    Real RhoH=-1.e20; pp.query("RhoH",RhoH);
    ostate.setVal(RhoH,box,sCompH,1);


    cout << "State Before: " << endl;
    Array<Real> probeVals(nSpec,0);
    for (int i=0; i<nSpec; i++) probeVals[i] = ostate(probe,sCompY+i);    

    if (do_sdc) {
      print_rstate(Temp, RhoH, Patm, probeVals, cd);
    }
    else {
      print_state(Temp, Patm, probeVals, cd);
    }
    
    FArrayBox funcCnt(box,1), tmp(box,nSpec+1);

    if (do_sdc) {
      std::cout << "Using SDC form of chemistry integrator..." << std::endl;
      ostate.setVal(0,sCompR);
      for (int i=0; i<nSpec; ++i)
      {
        ostate.plus(ostate,sCompY+i,sCompR,1);
      }

      for (int i=0; i<speciesNames.size(); ++i)
      {
        Real val = 0;
        std::string src_str = "c0_"+speciesNames[i];
        if (pp.countval(src_str.c_str())) {
          pp.query(src_str.c_str(),val);
          tmp.setVal(val,box,i,1);
        }
      }
      Real val = 0;
      std::string src_str = "c0_RhoH";
      pp.query(src_str.c_str(),val);
      tmp.setVal(val,box,speciesNames.size(),1);

      cd.solveTransient_sdc(nstate,nstate,nstate,ostate,ostate,ostate,tmp,funcCnt,
                            box,sCompY,sCompH,sCompT,dt,Patm,0);
      tmp.resize(nstate.box(),1);
      tmp.copy(ostate,sCompR,0,1);
      tmp.invert(1);
      for (int i=0; i<nSpec; ++i)
      {
        nstate.mult(tmp,0,sCompY+i,1);
      }
      nstate.mult(tmp,0,sCompH,1);     

      Real errMAX = 1.e-6;
      cd.getTGivenHY(nstate,nstate,nstate,box,sCompH,sCompY,sCompT,errMAX);

    } else {
      cd.solveTransient(nstate,nstate,ostate,ostate,funcCnt,
                        box,sCompY,sCompT,dt,Patm);
      
    }

    cout << "State After dt = " << dt << " sec: " << endl;
    for (int i=0; i<nSpec; i++) probeVals[i] = nstate(probe,sCompY+i);
    Temp = nstate(probe,sCompT);
    print_state(Temp, Patm, probeVals, cd);

    cout << " ... required " << funcCnt(probe,0) << " function evals" << endl;
    
    BoxLib::Finalize();
}

    
