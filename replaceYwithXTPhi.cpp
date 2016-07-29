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

using std::cout;
using std::cerr;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "FArrayBox.H"
#include "Utility.H"
#include "ChemDriver.H"
#include "WritePlotFile.H"

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  std::cout << '\t' << "[fuelname = <name of fuel species>]" << '\n';
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = BoxLib::Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    std::string fuelname="H2";     pp.query("fuelname",fuelname);

    Real Patm = 1.0;
    pp.query("Patm",Patm);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    ChemDriver cd;

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    bool output_fuel_only = false;
    pp.query("output_fuel_only",output_fuel_only);

    int do_heat_release(1);
    pp.query("do_heat_release",do_heat_release);
    
    int idSpIn = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "Y(" + cd.speciesNames()[0] + ")";
    const int nSpec = cd.numSpecies();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == spName) idSpIn = i;
    }
    if (ParallelDescriptor::IOProcessor() && idSpIn<0)
        cerr << "Cannot find required data in pltfile" << endl;

    const int idYlocal = 0;                   // Ys start here
    const int idTIn    = idYlocal+nSpec;      // Temp
    const int idRhoIn  = idTIn+1;             // Density
    const int idRhoHIn = idRhoIn+1;           // rho h
    const int idFCRIn  = idRhoHIn+1;          // FCR
    const int idHRIn   = idFCRIn+1;           // HR
    const int nCompIn  = idHRIn+do_heat_release; // Read in Species and Temp

    const int idSpOut  = 0;                   // Put Xs at start
    const int idT      = idSpOut+nSpec;       // Then Temp
    const int idPhi    = idT+1;               // Then Phi
    const int nCompOut = idPhi+1;             // Output Xs, Temp and Phi

    Array<std::string> outNames(nCompOut);
    Array<std::string> inNames(nCompIn);
    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nSpec; ++i) {
        destFillComps[i] = idYlocal + i;
        inNames[i] =  "Y(" + cd.speciesNames()[i] + ")";
        outNames[i] = "X(" + cd.speciesNames()[i] + ")";
    }
    for (int i=idYlocal+nSpec; i<nCompIn; i++) {
      destFillComps[i] = idYlocal + i;
    }

    inNames[idTIn]  = "temp";
    inNames[idRhoIn]  = "density";
    inNames[idRhoHIn]  = "rhoh";
    inNames[idFCRIn]  = fuelname+"_ConsumptionRate";
    if (do_heat_release==1)
	inNames[idHRIn]  = "HeatRelease";

    outNames[idT]   = "temp";
    outNames[idPhi] = "phi";
    
    PArray<MultiFab> outdata(Nlev,PArrayManage);
    PArray<MultiFab> indata(Nlev,PArrayManage);
    PArray<MultiFab> CData(Nlev,PArrayManage);
    PArray<MultiFab> CpData(Nlev,PArrayManage);
    const int nGrow = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        outdata.set(lev,new MultiFab(ba,nCompOut,nGrow));
        indata.set(lev,new MultiFab(ba,nCompIn,nGrow));
        CData.set(lev,new MultiFab(ba,nSpec,nGrow));
        CpData.set(lev,new MultiFab(ba,1,nGrow));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        amrData.FillVar(indata[lev],lev,inNames,destFillComps);
        for (int i=0; i<inNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        for (MFIter mfi(indata[lev]); mfi.isValid(); ++mfi)
        {
            const FArrayBox& Y = indata[lev][mfi];
            FArrayBox& X = outdata[lev][mfi];
            FArrayBox& C = CData[lev][mfi];
            FArrayBox& Cp = CpData[lev][mfi];
            const Box& box = mfi.validbox();

	    // Put Mole fractions in X
            cd.massFracToMoleFrac(X,Y,box,idYlocal,idSpOut);
	    // Put Temp in X
	    X.copy(Y,idTIn,idT,1);
	    // Put molar concentrations in C
            cd.massFracToMolarConc(C,Y,Y,Patm,box,idYlocal,idTIn,idSpOut); // Temp is in Y[idTIn]
	    FArrayBox C_C(box,1); if (fuelname!="H2") cd.getElementMoles(C_C,"C",C,box,0,0); else C_C.setVal(0.0);
            FArrayBox C_H(box,1); cd.getElementMoles(C_H,"H",C,box,0,0);
            FArrayBox C_O(box,1); cd.getElementMoles(C_O,"O",C,box,0,0);
#if 0
	    // This is the old way
	    if      ( fuelname=="H2"   ) {    X.copy(C_H,0,idPhi,1);    X.mult(0.5,idPhi,1);        X.divide(C_O,0,idPhi,1);}
	    else if ( fuelname=="CH4"  ) {    X.copy(C_C,0,idPhi,1);    X.mult(4.0,idPhi,1);        X.divide(C_O,0,idPhi,1);}
	    else if ( fuelname=="C3H8" ) {    X.copy(C_C,0,idPhi,1);    X.mult(10.0/3.0,idPhi,1);   X.divide(C_O,0,idPhi,1);}
#else
	    // This is the JBB way
	    C_H.mult(0.5);
	    C_C.mult(2.0);
	    C_C.plus(C_H);
	    C_C.divide(C_O);
	    X.copy(C_C,0,idPhi,1);
#endif
	    // Get Cp
	    cd.getCpmixGivenTY(Cp,Y,Y,box,idTIn,idYlocal,0);
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive finished for level " << lev << endl;
    }

    if (ParallelDescriptor::IOProcessor())
	cerr << fuelname << endl;	

  
    std::string nfile(getFileRoot(plotFileName) + "_XTPhi");

    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << nfile << endl;
    
    const bool verb = false;
    const AmrData& a = amrData;
    if (!output_fuel_only) {
      // APU
      WritePlotfile("NavierStokes-V1.1",outdata,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
		    a.DxLevel(),a.CoordSys(),nfile,outNames,verb);

    } else {
      // Just output X_fuel, Temp and Phi
      const int nCompOutFO = 11+do_heat_release;
      const int iFuel(cd.index(fuelname));
      Array<std::string> outNamesFO(nCompOutFO);
      PArray<MultiFab> outdataFO(Nlev,PArrayManage);
      outNamesFO[0] = "X(" + cd.speciesNames()[iFuel] + ")";
      outNamesFO[1] = "Y(" + cd.speciesNames()[iFuel] + ")";
      outNamesFO[2] = "C(" + cd.speciesNames()[iFuel] + ")";
      outNamesFO[3] = "temp";
      outNamesFO[4] = "phi";
      outNamesFO[5] = "rho";
      outNamesFO[6] = "rho.Y(" + cd.speciesNames()[iFuel] + ")";
      outNamesFO[7] = "rho.h";
      outNamesFO[8] = "FCR";
      if (do_heat_release==1)
	  outNamesFO[9] = "HR";
      outNamesFO[9+do_heat_release]= "Cp";
      outNamesFO[10+do_heat_release]= "CpT";
      
      for (int lev=0; lev<Nlev; ++lev)
	{
	  outdataFO.set(lev,new MultiFab(amrData.boxArray(lev),nCompOutFO,nGrow));
	  outdataFO[lev].copy(outdata[lev],   iFuel, 0,1); // X
	  outdataFO[lev].copy(indata[lev],    iFuel, 1,1); // Y
	  outdataFO[lev].copy(CData[lev],     iFuel, 2,1); // C
	  outdataFO[lev].copy(indata[lev],    idTIn, 3,1); // Temp
	  outdataFO[lev].copy(outdata[lev],   idPhi, 4,1); // Phi
	  outdataFO[lev].copy(indata[lev],  idRhoIn, 5,1); // rho
	  outdataFO[lev].copy(indata[lev],  idRhoIn, 6,1); // *rho*  Y
	  MultiFab::Multiply(outdataFO[lev],indata[lev],iFuel,6,1,nGrow); //  rho  *Y*
	  outdataFO[lev].copy(indata[lev], idRhoHIn, 7,1); // rho.h
	  outdataFO[lev].copy(indata[lev],  idFCRIn, 8,1); // FCR
	  if (do_heat_release==1)
	      outdataFO[lev].copy(indata[lev],   idHRIn, 9,1); // HR
	  outdataFO[lev].copy(CpData[lev],        0,9+do_heat_release,1); // Cp
	  outdataFO[lev].copy(CpData[lev],        0,10+do_heat_release,1); // *Cp*  T
	  MultiFab::Multiply(outdataFO[lev],indata[lev],idTIn,10+do_heat_release,1,nGrow); //  Cp  *T*
	}

      WritePlotfile("NavierStokes-V1.1",outdataFO,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
		    a.DxLevel(),a.CoordSys(),nfile,outNamesFO,verb);

    }

    BoxLib::Finalize();
    return 0;
}
