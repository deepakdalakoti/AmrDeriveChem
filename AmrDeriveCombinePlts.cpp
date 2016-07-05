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

#define Y_IN_PLOTFILE
#undef Y_IN_PLOTFILE

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infileL=<name> infileR=<name> outfile=<> [options] \n\tOptions:\n";
    exit(1);
}

void
clonePlotfile(AmrData& amrData,
              const PArray<MultiFab>& mfout,
              std::string& oFile,
              const Array<std::string>& names,
              bool verbose)
{
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(oFile,0755))
            BoxLib::CreateDirectoryFailed(oFile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    const int n_var = mfout[0].nComp();
    const int finestLevel = mfout.size() - 1;
    
    std::ofstream os;
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Opening file = " << oFileHeader << '\n';
    
        os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
        
        if (os.fail())
            BoxLib::FileOpenFailed(oFileHeader);
        //
        // Start writing plotfile.
        //
        os << amrData.PlotFileVersion() << '\n';
        os << n_var << '\n';
        for (int n = 0; n < n_var; n++) os << names[n] << '\n';
        os << BL_SPACEDIM << '\n';
        os << amrData.Time() << '\n';
        os << finestLevel << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) os << amrData.ProbLo()[i] << ' ';
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) os << amrData.ProbHi()[i] << ' ';
        os << '\n';
        for (int i = 0; i < finestLevel; i++) os << amrData.RefRatio()[i] << ' ';
        os << '\n';
        for (int i = 0; i <= finestLevel; i++) os << amrData.ProbDomain()[i] << ' ';
        os << '\n';
        for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
        os << '\n';
        for (int i = 0; i <= finestLevel; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << amrData.DxLevel()[i][k] << ' ';
            os << '\n';
        }
        os << amrData.CoordSys() << '\n';
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
        int nGrids = amrData.boxArray(iLevel).size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if (ParallelDescriptor::IOProcessor())
        {
            os << iLevel << ' ' << nGrids << ' ' << amrData.Time() << '\n';
            os << 0 << '\n';
            
            for (int i = 0; i < nGrids; ++i)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    os << amrData.GridLocLo()[iLevel][i][n]
                       << ' '
                       << amrData.GridLocHi()[iLevel][i][n]
                       << '\n';
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
        VisMF::Write(mfout[iLevel], PathName, VisMF::OneFilePerCPU);
    }
    
    if (ParallelDescriptor::IOProcessor()) {
      os.close();
    }
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

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    std::string infileL;
    pp.get("infileL",infileL);

    std::string infileR;
    pp.get("infileR",infileR);

    std::string outfile;
    pp.get("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Set up for reading left pltfile
    DataServices dataServicesL(infileL, fileType);

    if (!dataServicesL.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrDataL = dataServicesL.AmrDataRef();

    Array<int> compsL;
    if (int nc = pp.countval("compsL"))
    {
        compsL.resize(nc);
        pp.getarr("compsL",compsL,0,nc);
    }
    else
    {
        int sCompL = 0;
        pp.query("sCompL",sCompL);
        int nCompL = amrDataL.NComp();
        pp.query("nCompL",nCompL);
        BL_ASSERT(sCompL+nCompL <= amrDataL.NComp());
        compsL.resize(nCompL);
        for (int i=0; i<nCompL; ++i)
            compsL[i] = sCompL + i;
    }

    // Set up for reading right pltfile
    DataServices dataServicesR(infileR, fileType);

    if (!dataServicesR.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrDataR = dataServicesR.AmrDataRef();

    Array<int> compsR;
    if (int nc = pp.countval("compsR"))
    {
        compsR.resize(nc);
        pp.getarr("compsR",compsR,0,nc);
    }
    else
    {
        int sCompR = 0;
        pp.query("sCompR",sCompR);
        int nCompR = amrDataR.NComp();
        pp.query("nCompR",nCompR);
        BL_ASSERT(sCompR+nCompR <= amrDataR.NComp());
        compsR.resize(nCompR);
        for (int i=0; i<nCompR; ++i)
            compsR[i] = sCompR + i;
    }

    int Nlev = amrDataL.FinestLevel() + 1;
    BL_ASSERT(Nlev == amrDataR.FinestLevel() + 1);

    const int nComp = compsL.size() + compsR.size();
    PArray<MultiFab> fileData(Nlev);
    for (int lev=0; lev<Nlev; ++lev)
    {
        BL_ASSERT(amrDataL.boxArray(lev) == amrDataR.boxArray(lev));
        fileData.set(lev,new MultiFab(amrDataL.boxArray(lev),nComp,0));
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Full MultiFab allocated " << endl;

    for (int lev=0; lev<Nlev; ++lev)
    {
        for (int i=0; i<compsL.size(); ++i)
        {
            fileData[lev].copy(amrDataL.GetGrids(lev,compsL[i]),0,i,1);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "After GetGrids (L): " << amrDataL.PlotVarNames()[compsL[i]] << endl;
            //amrDataL.FlushGrids(compsL[i]);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "AmrData flushed (L): " << amrDataL.PlotVarNames()[compsL[i]] << endl;
        }
        for (int i=0; i<compsR.size(); ++i)
        {
            fileData[lev].copy(amrDataR.GetGrids(lev,compsR[i]),0,compsL.size()+i,1);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "After GetGrids (R): " << amrDataR.PlotVarNames()[compsR[i]] << endl;
            //amrDataR.FlushGrids(compsR[i]);
            if (ParallelDescriptor::IOProcessor())
                std::cerr << "AmrData flushed (R): " << amrDataR.PlotVarNames()[compsR[i]] << endl;
        }
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "File data loaded" << endl;

    Real progMin, progMax;
    for (int i=0; i<compsL.size(); ++i) {
        amrDataL.MinMax(amrDataL.ProbDomain()[Nlev-1], amrDataL.PlotVarNames()[compsL[i]], Nlev-1, progMin, progMax);
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << amrDataL.PlotVarNames()[compsL[i]] << " min/max: " << progMin << ", " << progMax << std::endl;
        }
    }
    for (int i=0; i<compsR.size(); ++i) {
        amrDataR.MinMax(amrDataR.ProbDomain()[Nlev-1], amrDataR.PlotVarNames()[compsR[i]], Nlev-1, progMin, progMax);
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << amrDataR.PlotVarNames()[compsR[i]] << " min/max: " << progMin << ", " << progMax << std::endl;
        }
    }

    Array<std::string> names(nComp);
    for (int i=0; i<compsL.size(); ++i)
        names[i] = amrDataL.PlotVarNames()[compsL[i]];
    for (int i=0; i<compsR.size(); ++i)
        names[compsL.size()+i] = amrDataR.PlotVarNames()[compsR[i]];
    
    bool verb = false;
    clonePlotfile(amrDataL,fileData,outfile,names,verb);

    BoxLib::Finalize();
    return 0;
}
