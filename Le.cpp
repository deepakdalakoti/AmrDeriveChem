#include <iomanip>

#include "REAL.H"
#include "Box.H"
#include "PArray.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"

using std::cout;
using std::endl;
using std::string;

int
main (int   argc,
      char* argv[])
{
    /*
      E.g., for c12 mixture at phi=0.7

      Le3d....ex SpecName=NC12H26 Y_NC12H26=.0448128430753 Y_O2=0.222479822384
    */
    BoxLib::Initialize(argc,argv);
    ParmParse pp;
    Real Patm=1; pp.query("Patm",Patm);
    Real Temp=298; pp.query("Temp",Temp);

    ChemDriver cd;
    int Nspec = cd.numSpecies();
    string SpecName = "NC12H26";
    pp.query("SpecName",SpecName);
    int idx = cd.index(SpecName);

    IntVect iv(D_DECL(0,0,0));
    Box bx(iv,iv);
    FArrayBox T(bx,1);     T.setVal(Temp);
    FArrayBox Y(bx,Nspec); Y.setVal(0);

    // Get mixture
    Real sum = 0;
    for (int i=0; i<Nspec; ++i) {
      string ppstr = "Y_" + cd.speciesNames()[i];
      pp.query(ppstr.c_str(),Y(iv,i));
      if (i == cd.index("N2") && Y(iv,i)!=0) {
	BoxLib::Abort("Do not set N2 mass fraction....I will do it");
      }
      sum += Y(iv,i);
    }
    Y(iv,cd.index("N2")) = 1 - sum;

    cout << "Temp: " << Temp << endl;
    cout << "Patm: " << Patm << endl;
    cout << "idx: " << idx << endl;
    for (int i=0; i<Nspec; ++i) {
      cout << "Y(" << cd.speciesNames()[i] << ") = " << Y(iv,i) << endl;      
    }

#ifdef LMC_SDC
    cout << "LMC_SDC" << endl;
    FArrayBox Rho(bx,1);
    cd.getRhoGivenPTY(Rho,Patm,T,Y,bx,0,0,0);
    Real rho = Rho(iv,0);
    cout << "Rho: " << rho << endl;
    for (int i=0; i<Nspec; ++i) {
      Y.mult(Rho,0,i,1);
    }
#endif

    bool do_temp = true;
    bool do_VelVisc = true;
    FArrayBox rhoD(bx,Nspec+2);
    cd.getMixAveragedRhoDiff(rhoD,Y,T,Patm,bx,0,0,0,do_temp,do_VelVisc);

    FArrayBox cpmix(bx,1);
    cd.getCpmixGivenTY(cpmix,T,Y,bx,0,0,0);

    Real lambda = rhoD(iv,Nspec);
    Real rD = rhoD(iv,idx);
    Real cp = cpmix(iv,0);
    Real Le = lambda/(rD*cp);

    cout << "lambda: " << lambda << endl;
    cout << "rD: " << rD << endl;
    cout << "cp: " << cp << endl;
    cout << "Le: " << Le << endl;

    BoxLib::Finalize();
    return 0;
}

