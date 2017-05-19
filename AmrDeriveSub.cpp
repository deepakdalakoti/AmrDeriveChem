//BL_COPYRIGHT_NOTICE

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
//using std::set_new_handler;
using std::cout;
using std::endl;

//#include <unistd.h>

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "WritePlotFile.H"
//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

const bool verbose_DEF   = false;

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    infile=inFileName" << '\n';
    cout << "   [outfile=outFileName <defaults to <inFileName>_section>]" << '\n';
    cout << "   [-sComp=N <defaults to 0, unless \"comps\" used]" << '\n';
    cout << "   [-nComp=N <defaults to all in <inFileName>, unless comps used]" << '\n';
    cout << "   [-comps=\"N1 N2 N3...\"]" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
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
//            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(data[iLevel], PathName, VisMF::OneFilePerCPU);
    }
    
    os.close();
}



static Array< Array<int> > contigLists(const Array<int> orig);

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);

    BoxLib::Initialize(argc,argv);    

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    Array<std::string> infiles(pp.countval("infiles"));

     // Read in pltfile and get amrData ref
     Array<int> comps;
    pp.getarr("infile",infiles);
    int nPlot = infiles.size();
   for (int iPlot= 0; iPlot < nPlot ; iPlot++) {
          std::string infile = infiles[iPlot];
    vector<std::string> pieces = BoxLib::Tokenize(infile,std::string("/"));
    std::string outfile = pieces[pieces.size()-1] + std::string("_section");
    pp.query("outfile",outfile);


    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = dataServices.AmrDataRef();
     if (int nc = pp.countval("comps"))
      {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
     }

    else

    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nComp = amrData.NComp();
        pp.query("nComp",nComp);
//        BL_ASSERT(sComp+nComp < amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }
    

    bool verbose = verbose_DEF;
    verbose = (pp.contains("verbose") ? true : false);
    if (verbose)
        AmrData::SetVerbose(true);

    Array< Array<int> > compsNEW = contigLists(comps);
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    Box subbox = amrData.ProbDomain()[finestLevel];
    Array<int> inBox;

    if (int nx=pp.countval("box"))
    {
        pp.getarr("box",inBox,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(inBox.size()==2*d);
        subbox=Box(IntVect(D_DECL(inBox[0],inBox[1],inBox[2])),
                   IntVect(D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
                   IndexType::TheCellType());
    }

    Array<Box> subboxes(finestLevel+1,subbox);
    for (int iLevel = finestLevel-1; iLevel>=0; --iLevel)
        subboxes[iLevel] = BoxLib::coarsen(subboxes[iLevel+1],amrData.RefRatio()[iLevel]);
    for (int iLevel = 1; iLevel<=finestLevel; ++iLevel)
        subboxes[iLevel] = BoxLib::refine(subboxes[iLevel-1],amrData.RefRatio()[iLevel-1]);

    PArray<MultiFab> data_sub(finestLevel+1);
    Array<std::string> names(comps.size());
    Array<std::string> subNames;
    Array<int> fillComps;
   
    Array<Real> plo(BL_SPACEDIM), phi(BL_SPACEDIM);
    Array<Box> psize(finestLevel+1);
    const IntVect ilo = subbox.smallEnd();
    const IntVect ihi = subbox.bigEnd();

   for (int i =0 ; i< BL_SPACEDIM; i++) {
       
       plo[i] = amrData.ProbLo()[i]+(ilo[i])*amrData.DxLevel()[finestLevel][i];
       phi[i] = amrData.ProbLo()[i]+(ihi[i]+1)*amrData.DxLevel()[finestLevel][i];
       
     }
        
    int actual_lev = -1;
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        // Build the BoxArray for the result
        BoxArray ba_sub = BoxLib::intersect(amrData.boxArray(iLevel),subboxes[iLevel]);

        if (ba_sub.size() > 0)
        {
             actual_lev = iLevel;
           data_sub.set(iLevel, new MultiFab(ba_sub,comps.size(),0,Fab_allocate));

            for (int i=0; i<comps.size(); ++i)
            {
                data_sub[iLevel].copy(amrData.GetGrids(iLevel,comps[i],subboxes[iLevel]),0,i,1);
                amrData.FlushGrids(comps[i]);                
                names[i] = amrData.PlotVarNames()[comps[i]];
                
                if (ParallelDescriptor::IOProcessor()) 
                    cout << "Filling " << names[i] << " on level " << iLevel << endl;
                
            }
        }
    }

    // Write out the subregion pltfile
//    WritePlotfile("NavierStokes-V1.1",data_sub,amrData.Time(),plo,phi,
//                  amrData.RefRatio(),subboxes,amrData.DxLevel(),amrData.CoordSys(),
//                  outfile,names,verbose);
   writePlotfile(data_sub,amrData.Time(),plo,phi,
                  amrData.RefRatio(),subboxes,amrData.DxLevel(),amrData.CoordSys(),
                  outfile,names,verbose);

   }//iPlot
    BoxLib::Finalize();
    return 0;
}

static Array< Array<int> > contigLists(const Array<int> orig)
{
    BL_ASSERT(orig.size() > 0);
    Array< Array<int> > res(1);
    int mySet = 0;
    res[mySet].resize(1);
    int myCnt = 0;
    res[mySet][myCnt] = orig[0];

    for (int i=1; i<orig.size(); ++i)
    {
        if (res[mySet][myCnt]+1 == orig[i])
        {
            res[mySet].resize(++myCnt + 1);
            res[mySet][myCnt] = orig[i];
        }
        else
        {
            myCnt = 0;
            res.resize(++mySet+1);
            res[mySet].resize(1);
            res[mySet][myCnt] = orig[i];
        }
    }
    return res;
}


