
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"

using std::cout;
using std::endl;
using std::cerr;

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    ParmParse pp;

    Real Patm=1; pp.query("Patm",Patm);
    Real Temp=298; pp.query("Temp",Temp);
    Real phi; pp.get("phi",phi);
    std::string fuelName; pp.get("fuelName",fuelName);

    ChemDriver cd;

    int nSpec = cd.numSpecies();

    // hc + a.O2 -> b.CO2 + c.H2O
    //for hc = CH3OCH3, a=3, b=2, c=3
    //for hc = CH4, a=2, b=1, c=2
    //for hc = H2, a=0.5, b=0, c=1

    Real a, b, c;
    if (fuelName=="CH4") {
        a = 2;
        b = 1;
        c = 2;
    }
    else if (fuelName == "H2") {
        a = 0.5;
        b = 0;
        c = 1;
    }
    else if (fuelName == "CH3OCH3") {
        a = 3;
        b = 4;
        c = 3;
    }
    else
        BoxLib::Abort("Do not yet know about this fuel");

    cout << "For " << fuelName << " at phi =" << phi << ":" << endl;
    cout << "  Temp = " << Temp << endl;
    cout << "  Patm = " << Patm << endl;

    Real XO2 = 1.0/(1.0 + phi/a  + 0.79/0.21);
    Real Xfu = phi * XO2 / a;
    Real XN2 = (0.79/0.21) * XO2;
    cout << "    X["<< fuelName << "] = " << Xfu << endl;
    cout << "     X[O2] = " << XO2 << endl;
    cout << "     X[N2] = " << XN2 << endl;

    std::cout << "Number of species: " << nSpec << std::endl;

    Array<Real> X(nSpec,0);
    X[cd.index(fuelName)] = Xfu;
    X[cd.index("O2")] = XO2;
    X[cd.index("N2")] = XN2;


    Array<Real> Y = cd.moleFracToMassFrac(X);
    cout << "    Y["<< fuelName << "] = " << Y[cd.index(fuelName)] << endl;
    cout << "     Y[O2] = " << Y[cd.index("O2")] << endl;
    cout << "     Y[N2] = " << Y[cd.index("N2")] << endl;

    IntVect iv(D_DECL(0,0,0));
    Box box(iv,iv);
    FArrayBox T(box,1);
    FArrayBox Yf(box,nSpec);

    T.setVal(Temp);
    Yf.setVal(0.0);
    Yf.setVal(Y[cd.index(fuelName)],cd.index(fuelName));
    Yf.setVal(Y[cd.index("O2")],cd.index("O2"));
    Yf.setVal(1 - Y[cd.index(fuelName)] - Y[cd.index("O2")],cd.index("N2"));

    FArrayBox rho(box,1);
    cd.getRhoGivenPTY(rho,Patm,T,Yf,box,0,0,0);
    cout << "density =" << rho(iv,0) << " kg/m^3" << endl;

    FArrayBox mu(box,1);
    cd.getMixShearVisc(mu,Yf,T,box,0,0,0);
    cout << "visc: " << mu(iv,0) << endl;

    FArrayBox rhoD(box,cd.numSpecies());
    cd.getMixAveragedRhoDiff(rhoD,Yf,T,Patm,box,0,0,0);
    cout << "D(fuel): " << rhoD(iv,cd.index(fuelName))/rho(iv,0) << endl;
    
    FArrayBox cpmix(box,1);
    cd.getCpmixGivenTY(cpmix,T,Yf,box,0,0,0);
    cout << "cpmix: " << cpmix(iv,0) << endl;

    cout << "1/(rY_fu): " << 1/(rho(iv,0)*Y[cd.index(fuelName)]) << endl;

    BoxLib::Finalize();
    return 0;
}



