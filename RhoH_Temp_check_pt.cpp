#include <iostream>
#include <cstdio>
using std::cout;
using std::endl;
using std::string;

#include "Utility.H"
#include "ParallelDescriptor.H"
#include "ChemDriver.H"
#include "ParmParse.H"
#include "VisMF.H"

const bool    verbose_DEF = true;
const string tranfile_DEF = "tran.asc.chem-H";

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    // Parse command line
    ParmParse pp;

    std::string tranfile=tranfile_DEF;   pp.query("tranfile",tranfile);
    ChemDriver cd(tranfile);

    const int nSpec = cd.numSpecies();
    Box box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    FArrayBox state(box,nSpec+2);
    int idT=nSpec;
    int idY=0;
    int idH=idT+1;
    state.setVal(0.999961490191e3,box,idT);
    state.setVal(0.538446978230e-3,idY+0);
    state.setVal(0.254524594123e-5,idY+1);
    state.setVal(0.875314236073e-4,idY+2);
    state.setVal(0.177067600194e+0,idY+3);
    state.setVal(0.857694946320e-4,idY+4);
    state.setVal(0.617993260292e-1,idY+5);
    state.setVal(0.671447280277e-4,idY+6);
    state.setVal(0.167183702821e-4,idY+7);
    state.setVal(0.760334917536e+0,idY+8);
    state.setVal(-0.241168225942e+09 * 1.e-4,idH);

    Real errMAX = 1.e-50; pp.query("errMAX",errMAX);

    cd.getTGivenHY(state,state,state,box,idH,idY,idT,errMAX);

    BoxLib::Finalize();
}

    
