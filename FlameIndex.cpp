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
#include "Geometry.H"
#include "FlameIndex_F.H"

#include <AppendToPlotFile.H>

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile plotFileNames=f1 f2 <where time(f1)<time(f2)> [options] \n\tOptions:\n";
    exit(1);
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

    const Real strt_time = ParallelDescriptor::second();
    Real strt_io, io_time = 0;

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    int verbose=0; pp.query("verbose",verbose);
    if (verbose>1)
        AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    strt_io = ParallelDescriptor::second();
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();
    io_time += ParallelDescriptor::second() - strt_io;
    ChemDriver cd; 
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    int nSpec = cd.numSpecies();
    const int nCompIn = nSpec+1+1;
    int indexO2,indexO, indexOH;
    Array<std::string> inVarAux(1);
    inVarAux[0] = "HeatRelease";
    Array<int> destAux(1);
    destAux[0] = 0;
    Array<std::string> inVarNames(nCompIn-1);
    for (int i=0; i<nSpec ; ++i) {
       inVarNames[i] = "Y("+cd.speciesNames()[i] + ")";
//       std::cout << cd.speciesNames()[i] <<std::endl; 
       if(cd.speciesNames()[i]=="O2") indexO2=i;
       if(cd.speciesNames()[i]=="O") indexO =i;
       if(cd.speciesNames()[i]=="OH") indexOH =i;
     }
   inVarNames[nSpec] = "HeatRelease";
//    inVarNames[idCst] = plotVarNames[idC];
    Array<int> destFillComps(nCompIn-1);
    for (int i=0; i<nCompIn-1; ++i)
        destFillComps[i] = i;
     
   
    Array<std::string> elements= cd.elementNames();
    int indC, indH, indO;
    for (int ind=0; ind< elements.size(); ind++)
    {
         if(elements[ind]=="c" || elements[ind] =="C")
              indC = ind;
         if(elements[ind]=="h" || elements[ind] =="H")
              indH = ind;
         if(elements[ind]=="o" || elements[ind] =="O")
              indO = ind;


  }
    const int nCompOut = 1 ; // 1 component stores the ||gradT||
    PArray<MultiFab> Specs(Nlev,PArrayManage);
    PArray<MultiFab> Out(Nlev,PArrayManage);
    PArray<MultiFab> Aux(Nlev,PArrayManage);

    PArray<Geometry> geoms(Nlev,PArrayManage);
    const int nGrow = 1;
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "pltfile has Nlev: " << Nlev << std::endl;
    }
     Real Max = 0;
    for (int lev=0; lev<Nlev; ++lev)
    {
        strt_io = ParallelDescriptor::second();
        const BoxArray ba = amrData.boxArray(lev);
        Out.set(lev,new MultiFab(ba,1,0));
        Aux.set(lev,new MultiFab(ba,1,0)); 
            
        Specs.set(lev, new MultiFab(ba,nCompIn,nGrow));
        const Array<Real>& delta = amrData.DxLevel()[lev];

        // Get input state data onto intersection ba
        const int myNComp = nCompIn; // gonna need this for fortran calls

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;
 //       Specs[lev].setVal(0);
        Out[lev].setVal(0);
        Specs[lev].setBndry(0.0); // Initialize grow cells to 0....here don't care to get grads right on c-f
        Specs[lev].setVal(0.0);
        Aux[lev].setVal(0.0);
        ParallelDescriptor::Barrier();
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Doing FillVar at lev: " << lev << std::endl;
        }
        amrData.FillVar(Specs[lev],lev,inVarNames,destFillComps);
        amrData.FillVar(Aux[lev],lev,inVarAux,destAux);
       for (int k=0; k<nSpec ; ++k) {
         int nC = cd.numberOfElementXinSpeciesY(elements[indC],cd.speciesNames()[k]);
         int nH = cd.numberOfElementXinSpeciesY(elements[indH],cd.speciesNames()[k]);
         int nO = cd.numberOfElementXinSpeciesY(elements[indO],cd.speciesNames()[k]);

           if((nC ==12 && nH ==26) || (nC==1 && nO ==1 && nH == 0 )|| (nC==12 && nO ==2 && nH == 25)) {
            if(ParallelDescriptor::IOProcessor())
                std::cout << " Species added " << cd.speciesNames()[k] << std::endl;
            MultiFab::Add(Specs[lev],Specs[lev],k,nSpec+1,1,1);
                }
            }
        MultiFab::Copy(Specs[lev],Specs[lev],indexO2,nSpec,1,1);
//        MultiFab::Add(Specs[lev],Specs[lev],indexO,nSpec,1,1);
//        MultiFab::Add(Specs[lev],Specs[lev],indexOH,nSpec,1,1);


        ParallelDescriptor::Barrier();

        io_time += ParallelDescriptor::second() - strt_io;
        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        const bool do_corners = false;
        geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));

        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(Specs[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = Specs[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                          s.dataPtr(0),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                          &myNComp);
        }

        // Fix up fine-fine and periodic for idCst
        Specs[lev].FillBoundary(0,nCompIn);
//        geoms[lev].FillPeriodicBoundary(Specs[lev],0,nCompIn,do_corners);

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Boundary cells filled at lev: " << lev << std::endl;
        }
        // Compute mag of progress gradient.  Result in state, comp=idGr
        FArrayBox nWork;
        for (MFIter mfi(Specs[lev]); mfi.isValid(); ++mfi)
        {
            
            FArrayBox& s = Specs[lev][mfi];
            const Box& box = mfi.validbox();
            FArrayBox& OutBox = Out[lev][mfi];    
            FORT_GRAD(box.loVect(),box.hiVect(),
                      s.dataPtr(nSpec+1),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                      s.dataPtr(nSpec), ARLIM(s.loVect()),ARLIM(s.hiVect()),
                      delta.dataPtr(),OutBox.dataPtr(),ARLIM(OutBox.loVect()),
                      ARLIM(OutBox.hiVect()),Aux[lev][mfi].dataPtr(0));
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Gradient vector has been computed on level " << lev << endl;

        ParallelDescriptor::Barrier();

        Max  = std::max(Out[lev].max(0,0),Max);
    }
    
//     for (int lev =0 ; lev < Nlev ; lev ++)   // normalize by maxima
 ///         Out[lev].mult(1/Max);

    /* adding one component in the nCompOut to get the ||grad progress ||, and this might be useful to compute the thermal thickness. */
    Array<std::string> nnames(nCompOut);
    nnames[0]= "FlameIndex_4";
    bool verb=false;
    bool appendPlotFile=false; pp.query("append",appendPlotFile);

    strt_io = ParallelDescriptor::second();
       std::string outfile(getFileRoot(plotFileName) + "_FlameIndex"); pp.query("outfile",outfile);

        if (ParallelDescriptor::IOProcessor())
            cout << "Appending new data to " << plotFileName << endl;

        const AmrData& a = amrData;
        ParallelDescriptor::Barrier();
//        Array<std::string> nnames(1);
        nnames[0] = "FlameIndex_4";
       if (!appendPlotFile) {
        writePlotfile(Out,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),

                      a.DxLevel(),a.CoordSys(),outfile,nnames,verb);
             }
        else {
        ParallelDescriptor::Barrier();
      
       appendToPlotFile(amrData, Out, plotFileName, nnames,"FlameIndex_4",0);
       }

    io_time += ParallelDescriptor::second() - strt_io;

    //long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
    //long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
    
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
/*
    ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nFAB byte spread across MPI nodes: ["
                  << min_fab_bytes
                  << " ... "
                  << max_fab_bytes
                  << "]\n";
    }
*/
    const Real end_time = ParallelDescriptor::second();
    const Real run_time = end_time - strt_time - io_time;

    Real run_time_max, run_time_min; run_time_max = run_time_min = run_time;
    Real io_time_max, io_time_min; io_time_max = io_time_min = io_time;

    ParallelDescriptor::ReduceRealMax(run_time_max,IOProc);
    ParallelDescriptor::ReduceRealMax(io_time_max,IOProc);
    ParallelDescriptor::ReduceRealMin(run_time_min,IOProc);
    ParallelDescriptor::ReduceRealMin(io_time_min,IOProc);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Max Compute time: " << run_time_max << '\n';
        std::cout << "Min Compute time: " << run_time_min << '\n';
        std::cout << "Max I/O time: " << io_time_max << '\n';
        std::cout << "Min I/O time: " << io_time_min << '\n';
    }


    BoxLib::Finalize();
    return 0;
}
