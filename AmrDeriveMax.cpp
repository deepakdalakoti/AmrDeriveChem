//BL_COPYRIGHT_NOTICE

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
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

    if (int nx=pp.countval("box"))  //in finest level indexes where finest level is the level provided
    {
        pp.getarr("box",inBox,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(inBox.size()==2*d);
        subbox=Box(IntVect(D_DECL(inBox[0],inBox[1],inBox[2])),
                   IntVect(D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
                   IndexType::TheCellType());
    }


    PArray<MultiFab> data_sub(1);
    Array<std::string> names(comps.size());
    Array<std::string> subNames;
    Array<int> fillComps(comps.size());
    int indZ;
    for (int i =0 ; i< comps.size() ; i++)
    {   names[i] = amrData.PlotVarNames()[comps[i]];
       if (names[i]=="Z") indZ  =i;
       fillComps[i]=i;
     }

     
    Array<Real> plo(BL_SPACEDIM), phi(BL_SPACEDIM);
    Array<Box> psize(finestLevel+1);
    const IntVect ilo = subbox.smallEnd();
    const IntVect ihi = subbox.bigEnd();

   for (int i =0 ; i< BL_SPACEDIM; i++) {
       
       plo[i] = amrData.ProbLo()[i]+(ilo[i])*amrData.DxLevel()[finestLevel][i];
       phi[i] = amrData.ProbLo()[i]+(ihi[i]+1)*amrData.DxLevel()[finestLevel][i];
       
     }
   Array<Real> dx = amrData.DxLevel()[finestLevel];  
   Box domain = amrData.ProbDomain()[finestLevel];
   Array<int> lo(3);
   Array<int> hi(3);
   for (int i=0; i < 3 ; i ++ )
   {
        lo[i] = domain.smallEnd()[i];
        hi[i] = domain.bigEnd()[i];
  }
//  std::cout << ilo[0] <<" " <<  ihi[0] << std:: endl;
//    BoxArray ba = amrData.boxArray(finestLevel);
    BoxArray ba = amrData.boxArray(0);
    int Nlev = amrData.FinestLevel()+1;
    for (int i=0 ; i < Nlev-1 ; i++)
       ba.refine(2);

    BoxArray ba_i = BoxLib::intersect(ba,subbox);
    Array<Box> outB(1);
       outB[0] = subbox;
    data_sub.set(0,new MultiFab(ba_i,comps.size(),0,Fab_allocate));
    amrData.FillVar(data_sub[0],finestLevel,names,fillComps);
   Array<Real> MaxT(hi[2]+1);
   for (MFIter mfi(data_sub[0]); mfi.isValid() ; ++mfi )  {
     int proc = ParallelDescriptor::MyProc();

      const Box& vbox = mfi.validbox();
      int zlo = vbox.smallEnd()[2];
      int zhi = vbox.bigEnd()[2];
      for (int k = zlo; k < zhi ; ++k) {
//      Box maxbox = Box(IntVect(D_DECL(lo[0],lo[1],k)),IntVect(D_DECL(hi[0],hi[1],k)));
//      if( k > 2558 )
//        continue;
      MaxT[k] = std::max(MaxT[k], data_sub[0][mfi].max(0));
//      std::cout << data_sub[0][mfi].max(0) << std::endl;
//      MaxT[k] = 1.0;

      ParallelDescriptor::Bcast(&MaxT[k],1,proc);       
   }
     

//                      samplefile.close(); move this below
                   }



      if (ParallelDescriptor::IOProcessor()) {

	     
 	  FILE *file;
	  char filename[512];
          sprintf(filename,"%s", "MaxTemp.dat") ;
	  file = fopen(filename,"w");
          std::cout << "      filename = " << filename << std::endl;
	      for (int i=0; i<hi[2]; i++) {
		  fprintf(file,"%e ",MaxT[i]);
                  fprintf(file,"\n");

              }
	  fclose(file);
      
}

    // Write out the subregion pltfile

/*
  Array<int> Ref(1);
  Ref[0]= amrData.RefRatio()[0];
  Array<Array<Real> > dx1(1);
  dx1[0].resize(3);
   dx1[0] = amrData.DxLevel()[finestLevel];
   writePlotfile(data_sub,amrData.Time(),plo,phi,
                  Ref,outB,dx1,amrData.CoordSys(),
                  outfile,names,verbose);
*/
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


