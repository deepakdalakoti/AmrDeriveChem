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

#include "Geometry.H"
#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include "budget2_F.H"

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
    std::cerr << "\t TransportFile=[string]\n";
    std::cerr << "\t finestLevel=[int]\n";
    std::cerr << "\t verbose=[int]\n";
    exit(1);
}

void
writePlotfile(const PArray<MultiFab>&    data,
              Real                       time,
              const Array<Real>&         probLo,
              const Array<Real>&         probHi,
              const Array<int>&          refRatio,
              const Array<Box>&          probDomain,
              const Array<Array<Real> >& dxLevel,
              int                        coordSys,
              std::string&               oFile,
              const Array<std::string>&  names,
              bool                       verbose)
{
    // This is the version of plotfile that will be written
    std::string plotFileVersion = "NavierStokes-V1.1";

    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(oFile,0755))
            BoxLib::CreateDirectoryFailed(oFile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    std::ofstream os;
    const int finestLevel = data.size() - 1;

    if (ParallelDescriptor::IOProcessor())
    {

        std::string oFileHeader(oFile);
        oFileHeader += "/Header";
        
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        
        //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        
        if (verbose)
            std::cout << "Opening file = " << oFileHeader << '\n';
        
        os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
        
        if (os.fail())
            BoxLib::FileOpenFailed(oFileHeader);
        //
        // Start writing plotfile.
        //
        os << plotFileVersion << '\n';
        int n_var = data[0].nComp();
        os << n_var << '\n';
        for (int n = 0; n < n_var; n++) os << names[n] << '\n';
        os << BL_SPACEDIM << '\n';
        os << time << '\n';
        os << finestLevel << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
        os << '\n';
        for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
        os << '\n';
        for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
        os << '\n';
        for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
        os << '\n';
        for (int i = 0; i <= finestLevel; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << dxLevel[i][k] << ' ';
            os << '\n';
        }
        os << coordSys << '\n';
        os << "0\n"; // The bndry data width.
    }

    //
    // Write out level by level.
    //
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Write state data.
        //
        const BoxArray& ba = data[iLevel].boxArray();
        int nGrids = ba.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if (ParallelDescriptor::IOProcessor())
        {
            os << iLevel << ' ' << nGrids << ' ' << time << '\n';
            os << 0 << '\n';
            
            for (int i = 0; i < nGrids; ++i)
            {
                const Box& b = ba[i];
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    Real glo = b.smallEnd()[n]*dxLevel[iLevel][n];
                    Real ghi = (b.bigEnd()[n]+1)*dxLevel[iLevel][n];
                    os << glo << ' ' << ghi << '\n';
                }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            //
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
            if (!BoxLib::UtilCreateDirectory(Level, 0755))
                BoxLib::CreateDirectoryFailed(Level);
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        static const std::string MultiFabBaseName("/MultiFab");
        
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += MultiFabBaseName;
        
        if (ParallelDescriptor::IOProcessor())
        {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(data[iLevel], PathName, VisMF::OneFilePerCPU);
    }
    
    os.close();
}

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
            
            if (!BoxLib::UtilCreateDirectory(Level, 0755))
                BoxLib::CreateDirectoryFailed(Level);
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

            std::string oldMFname;
            is >> oldMFname;
            os << oldMFname << '\n';
            os << RelativePathName << '\n';
        }
        VisMF::Write(mfout[iLevel], PathName, VisMF::OneFilePerCPU);
    }
    
    os.close();
    is.close();
}

vector<std::string>
tokenize (const std::string& instr, const std::string& separators)
{
    vector<char*> ptr;
    //
    // Make copy of line that we can modify.
    //
    char* line = new char[instr.size()+1];

    (void) strcpy(line, instr.c_str());

    char* token = 0;

    if (!((token = strtok(line, separators.c_str())) == 0))
    {
        ptr.push_back(token);
        while (!((token = strtok(0, separators.c_str())) == 0))
            ptr.push_back(token);
    }

    vector<std::string> tokens(ptr.size());

    for (int i = 1; i < ptr.size(); i++)
    {
        char* p = ptr[i];

        while (strchr(separators.c_str(), *(p-1)) != 0)
            *--p = 0;
    }

    for (int i = 0; i < ptr.size(); i++)
        tokens[i] = ptr[i];

    delete line;

    return tokens;
}

std::string
getFileRoot(const std::string& infile)
{
    vector<std::string> tokens = tokenize(infile,std::string("/"));
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

    int verbose=0; pp.query("verbose",verbose);
    if (verbose>1)
        AmrData::SetVerbose(true);

    ChemDriver cd;

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();


    int idXin = -1;
    int idTin = -1;
    int idZ =-1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "Y(" + cd.speciesNames()[0] + ")";
    const int nSpec = cd.numSpecies();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == spName) idXin = i;
        if (plotVarNames[i] == "temp") idTin = i;
        if (plotVarNames[i] == "Z") idZ =i;
    }
    if (ParallelDescriptor::IOProcessor() && (idXin<0 || idTin<0) )
        cerr << "Cannot find required data in pltfile" << endl;

    Array<std::string> var2(1);
    var2[0] = "density";
    Array<int> dest2(1);
    dest2[0]=0;
    if(ParallelDescriptor::IOProcessor()) std::cout << "checkpoint 1 " << std::endl;
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    Array<std::string> whichSpec(pp.countval("whichSpec"));
    pp.getarr("whichSpec",whichSpec);
 //   std::cout << "which Spec "<< whichSpec[0] << whichSpec.size() << std::endl;
    int Nlev = finestLevel + 1;
    Array<int> NwhichSpec(whichSpec.size());
    for (int i=0; i<whichSpec.size(); ++i) {
        NwhichSpec[i] = cd.index(whichSpec[i]);
    }
    if(ParallelDescriptor::IOProcessor()) std::cout << "checkpoint 2 " << std::endl;

    const int idXst = 0;
    const int idTst = nSpec;
    const int nCompIn = idTst + 1+BL_SPACEDIM+1;
    const int si = whichSpec.size();
 if(ParallelDescriptor::IOProcessor()) std::cout << "checkpoint 3 " << std::endl;

    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;
    Array<std::string> inVarNames(nCompIn);
    for (int i=0; i<nSpec; ++i){
        inVarNames[idXst+i] = plotVarNames[idXin+i];
 if(ParallelDescriptor::IOProcessor()) std::cout << inVarNames[idXst+i] <<std::endl; }
    inVarNames[idTst] = plotVarNames[idTin];
    inVarNames[idTst+1] ="density"; 
    inVarNames[idTst+2] = "x_velocity";
    inVarNames[idTst+3] = "y_velocity";
#if BL_SPACEDIM==3
    inVarNames[idTst+4] = "z_velocity";
#endif

    PArray<MultiFab> outState(Nlev,PArrayManage);
    PArray<Geometry> geoms(Nlev,PArrayManage);
    PArray<MultiFab> RHOD(Nlev, PArrayManage); 
    const int nGrow = 1;

    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        MultiFab inState(ba,nCompIn,nGrow);
        MultiFab rho(ba,1,0); 
        MultiFab MixF(ba,1,0);
        outState.set(lev,new MultiFab(ba,8,0)); //1 for mixture frac
        RHOD.set(lev, new MultiFab(ba,nSpec,1));
        Array<Real> dx(BL_SPACEDIM);
        for (int i=0; i<BL_SPACEDIM; ++i)
            dx[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        amrData.FillVar(inState,lev,inVarNames,destFillComps);
        amrData.FillVar(rho,lev,var2,dest2);
        amrData.FillVar(MixF,lev,plotVarNames[idZ],0);
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        const bool do_corners = false;
        geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));

        // Extrap grow cells, as a first guess
        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(inState); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = inState[mfi];
            const Box& box = mfi.validbox();
            FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                          s.dataPtr(),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                          &nCompIn);
        }

        // Fix up fine-fine and periodic
        inState.FillBoundary(0,nCompIn);
        geoms[lev].FillPeriodicBoundary(inState,0,nCompIn,do_corners);

        FArrayBox RhoD,Reac,rhoY,Y,DTerms,rhoH;
        for (MFIter mfi(inState); mfi.isValid(); ++mfi)
        {
            const Box& vbox = mfi.validbox();
            const Box gbox = Box(vbox).grow(nGrow);
            FArrayBox AdvT;          
            
            RhoD.resize(gbox,nSpec); // last comp will hold lambda
            Reac.resize(vbox,nSpec);
            rhoY.resize(vbox,nSpec);
            rhoH.resize(vbox,1);
            AdvT.resize(vbox,3);
            Y.resize(vbox,nSpec);
            DTerms.resize(vbox,3*nSpec);
            RhoD.setVal(0.0);
            rhoY.setVal(0.0);
            Reac.setVal(0.0);
            DTerms.setVal(0.0);
            rhoH.setVal(0.0);
            AdvT.setVal(0.0);
            const FArrayBox& mixf = MixF[mfi];
            const FArrayBox& Y = inState[mfi];
            const FArrayBox& T = inState[mfi];
            FArrayBox& OutB = outState[lev][mfi];
            OutB.setVal(0.0);
            rhoY.copy(Y,0,0,nSpec);
            for (int i=0 ; i<nSpec; i++) {
              rhoY.mult(T,nSpec+1,i,1);
              }
            const Real Patm = 60.0; 
            const int do_temp = 0;
            const int do_VelVisc = 0;
            FORT_MIXAVG_RHODIFF_TEMP(gbox.loVect(),gbox.hiVect(),
                                     RhoD.dataPtr(),ARLIM(RhoD.loVect()),ARLIM(RhoD.hiVect()),
                                     T.dataPtr(idTst),ARLIM(T.loVect()),ARLIM(T.hiVect()),
                                     rhoY.dataPtr(),ARLIM(rhoY.loVect()),ARLIM(rhoY.hiVect()),
                                     
                                     &Patm, &do_temp, &do_VelVisc);
          
             
              FArrayBox& RHODBox = RHOD[lev][mfi];
              RHODBox.copy(RhoD,0,0,nSpec);
             

            FORT_COMPDIFFTERMS(vbox.loVect(),vbox.hiVect(),
                               Y.dataPtr(),         ARLIM(Y.loVect()),      ARLIM(Y.hiVect()),
                            RhoD.dataPtr(),      ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                             T.dataPtr(idTst+1), ARLIM(T.loVect()), ARLIM(T.hiVect()),
                               T.dataPtr(idTst+2), ARLIM(T.loVect()),ARLIM(T.hiVect()),
                               DTerms.dataPtr(),    ARLIM(DTerms.loVect()) ,ARLIM(DTerms.hiVect()),
                               dx.dataPtr(), &nSpec, &si,NwhichSpec.dataPtr(), AdvT.dataPtr());

          
            OutB.copy(DTerms,NwhichSpec[0],0,1);
            OutB.copy(DTerms,NwhichSpec[0]+nSpec,1,1);
            OutB.copy(DTerms,NwhichSpec[0]+nSpec+nSpec,2,1);
            cd.getHmixGivenTY(rhoH,T,Y,vbox,nSpec,0,0);
            rhoY.copy(Y,0,0,nSpec);
            rhoH.mult(Y,nSpec+1,0,1);
            for(int i=0; i<nSpec; i++) {
               rhoY.mult(Y,nSpec+1,i,1); 
                  }
//            OutB.copy(AdvT,0,4,1);
//            OutB.copy(AdvT,1,5,1);
//            OutB.copy(AdvT,2,6,1);

//            OutB.copy(mixf,0,7,1);

//      std::cout << "chekpoint 2" << std::endl;
          
          cd.reactionRateRhoY(Reac,rhoY,rhoH,T,Patm,vbox,0,0,nSpec,0);
//          cd.reactionRateY(Reac,Y,T,Patm,vbox,0,idTst,0);
//          Reac.mult(rho[mfi],0,NwhichSpec[0],1);

       for (int k=0 ; k<whichSpec.size() ; ++k) {
           OutB.copy(Reac,NwhichSpec[k],3,1);
        } 
   }       

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive completed on level " << lev << endl;
   
    }

    std::string outfile(getFileRoot(plotFileName) + "_budget"); pp.query("outfile",outfile);
    
    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << outfile << endl;
    
    Array<string> nnames((2*BL_SPACEDIM+1)*whichSpec.size()+1);
  for (int j=0 ; j<whichSpec.size(); ++j) {
   
        
          nnames[j*whichSpec.size()+0]= "Diff_X_" + whichSpec[j] ;
          nnames[j*whichSpec.size()+1]= "Diff_Y_" + whichSpec[j] ;
          nnames[j*whichSpec.size()+2] = "Diff_Z" + whichSpec[j];
          nnames[j*whichSpec.size()+3] = "Reaction_" +  whichSpec[j];
    }
          nnames[4] = "ADV_X";
          nnames[5] = "ADV_Y";
          nnames[6] = "ADV_Z";
          nnames[7] = "Z";


    bool verb=false;
    const AmrData& a = amrData;
    writePlotfile(outState,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),outfile,nnames,verb);
   

    
    BoxLib::Finalize();
    return 0;
}
