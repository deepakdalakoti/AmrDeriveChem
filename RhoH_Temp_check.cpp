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

    int vin=(verbose_DEF?1:0); pp.query("verbose",vin);
    bool verbose = (vin==1 ? true : false);
    
    std::string tranfile=tranfile_DEF;   pp.query("tranfile",tranfile);

    std::string fineFile="hackLev3State.mfab"; pp.query("fineFile",fineFile);
    std::string crseFile="hackLev2State.mfab"; pp.query("crseFile",crseFile);

    ChemDriver cd(tranfile);
    if (verbose)
        cd.set_verbose_vode();

    const Array<std::string>& names = cd.speciesNames();
    const int nSpec = cd.numSpecies();

    MultiFab mfC, mfF;
    VisMF::Read(mfF,fineFile);
    VisMF::Read(mfC,crseFile);

    const Box& bc = mfC.boxArray()[0];
    const Box& bf = mfF.boxArray()[0];

    int rat = 2;

    Box bcomC = Box(Box(bf).coarsen(rat) & bc);
    Box bcomF = Box(bcomC).refine(rat);

    int cnt = 0; 
    int x_velocity = cnt++;
    int y_velocity = cnt++;
    int z_velocity = cnt++;
    int Density = cnt++;
    int FirstSpec = cnt;
    cnt += nSpec;
    int RhoH = cnt++;
    int Trac = cnt++;
    int Temp = cnt++;
    int RhoRT = cnt++;
    int nComp = cnt;

    FArrayBox fabC(bcomC,nComp); fabC.copy(mfC[0],0,0,nComp);
    FArrayBox fabF(bcomF,nComp); fabF.copy(mfF[0],0,0,nComp);

    IntVect ivC(618,633,226);
    Box boxC(ivC,ivC);
    Box boxF = Box(boxC).refine(2);
    cout << "boxC: " << boxC << endl;
    cout << "boxF: " << boxF << endl;


    IntVect ivMinC = fabC.minIndex(bcomC,Temp);
    cout << "ivMinC: " << ivMinC << endl;
    cout << "min val: " << fabC(ivMinC,Temp) << endl;

    IntVect ivMinF = fabF.minIndex(bcomF,Temp);
    cout << "ivMinF: " << ivMinF << endl;
    cout << "min val: " << fabF(ivMinF,Temp) << endl;

    Real avg=0;
    cnt = 0;
    for (IntVect iv1=boxF.smallEnd(); iv1<=boxF.bigEnd(); boxF.next(iv1))
    {
        avg += fabF(iv1,Temp);
        cnt++;
    }
    avg *= 1./cnt;
    cout << " ... avg fine Temp: " << avg << endl;

    FArrayBox fabYHf(boxF,nSpec+1);
    fabYHf.copy(fabF,FirstSpec,0,nSpec);
    fabYHf.copy(fabF,RhoH,nSpec,1);
    FArrayBox newTf(boxF,1);
    for (int i=0; i<nSpec; ++i) {
        fabYHf.divide(fabF,Density,i,1);
    }
    fabYHf.divide(fabF,Density,nSpec,1);
    newTf.copy(fabF,Temp,0,1);
    cd.getTGivenHY(newTf,fabYHf,fabYHf,boxF,nSpec,0,0);

    cout << "old,new T: " << endl;
    for (IntVect iv1=boxF.smallEnd(); iv1<=boxF.bigEnd(); boxF.next(iv1))
    {
        cout << iv1 << " " << fabF(iv1,Temp) << " " << newTf(iv1,0) << endl;
    }
    BoxLib::Finalize();
}

    
