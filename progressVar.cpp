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
#include "progressVar_F.H"

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
              
              std::cout << "putback" << oldMFname[0] << std::endl;
             break;    }  
             test = oldMFname;
            std::cout << test << std::endl;     
            
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

    std::ifstream infile;
     Real x, y;
    Array< Array<Real> > TempZ(2);
    infile.open("Temp_Z_sorted.dat");
    if(!infile)
      {
           cerr << " Unable to open file "<<  cerr << endl ;

       }
    while(infile >> x >> y) {
         TempZ[0].push_back(x);
         TempZ[1].push_back(y);
         }
    
    int size = TempZ[0].size(); 
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
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    
    const int nSpec = cd.numSpecies();
    
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    const int nCompIn =2; // Z and temp
    
    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;
    Array<std::string> inVarNames(nCompIn);
    inVarNames[0] = "temp";
    inVarNames[1] = "Z";
    PArray<MultiFab> outState(Nlev,PArrayManage);
    
      for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        MultiFab inState(ba,nCompIn,0);
        outState.set(lev,new MultiFab(ba,1,0));
         Array<Real> dx(BL_SPACEDIM);
        for (int i=0; i<BL_SPACEDIM; ++i)
            dx[i] = amrData.ProbSize()[i] / amrData.ProbDomain()[lev].length(i);

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        amrData.FillVar(inState,lev,inVarNames,destFillComps);
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        
        for (MFIter mfi(inState); mfi.isValid(); ++mfi)
        {
           const Box& vbox = mfi.validbox();
          
            
            const FArrayBox& T = inState[mfi]; //first comp is temp second is Z
            FArrayBox& out  =    outState[lev][mfi];

            FORT_COMPDIFFTERMS(vbox.loVect(),vbox.hiVect(),
                               T.dataPtr(),         ARLIM(T.loVect()),      ARLIM(T.hiVect()),
                              
                               out.dataPtr(),TempZ[0].dataPtr(),TempZ[1].dataPtr(),&size);

   }       

      if (ParallelDescriptor::IOProcessor())
          cerr << "Derive completed on level " << lev << endl;

 }
  std::string outfile; pp.query("outfile",outfile);
  Array<std::string> nnames(1);
  nnames[0] = "ProgressVar";
 bool verb=false;
  const AmrData& a = amrData;
   bool append = true ; 
 
  if (ParallelDescriptor::IOProcessor())
      cout << "Writing new data to " << outfile << endl;
  

    pp.query("append", append);
   if(append) {
   appendToPlotFile(amrData, outState, plotFileName, nnames,"Damkoh",0);
    }
   else
   {
  writePlotfile(outState,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
               a.DxLevel(),a.CoordSys(),outfile,nnames,verb);
   }
  
    BoxLib::Finalize();
    return 0;
}
