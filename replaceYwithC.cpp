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

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
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

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

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

    int idYin = -1;
    int idTin = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName = "Y(" + cd.speciesNames()[0] + ")";
    const std::string TName = "temp";
    const int nSpec = cd.numSpecies();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == spName) idYin = i;
        if (plotVarNames[i] == TName) idTin = i;
    }
    if (ParallelDescriptor::IOProcessor() && (idYin<0 || idTin<0)) {
        BoxLib::Abort("Cannot find required data in pltfile");
    }

    const int idCout = 0;
    const int idTout = nSpec;
    const int nCompIn = nSpec+1;
    const int nCompOut = nSpec + 1;
    Real Patm = 1.0; pp.query ("Patm", Patm);
    
    Array<std::string> outNames(nCompOut);
    Array<std::string> inNames(nCompIn);
    Array<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Ys start here
    const int idTlocal = nSpec; // T starts here
    for (int i=0; i<nSpec; ++i)
    {
        destFillComps[i] = idYlocal + i;
        inNames[i] =  "Y(" + cd.speciesNames()[i] + ")";
        outNames[i] = "C(" + cd.speciesNames()[i] + ")";
    }
    destFillComps[idTlocal] = idTlocal;
    inNames[idTlocal] = TName;
    outNames[idTout] = TName;
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "To read:" << std::endl;
        for (int i=0; i<nCompIn; ++i) {
	  std::cout << inNames[i] << " into comp: " << destFillComps[i] << std::endl;
	}
        std::cout << "To write:" << std::endl;
        for (int i=0; i<nCompOut; ++i) {
	  std::cout << outNames[i] << " into comp: " << destFillComps[i] << std::endl;
	}
	std::cout << "T from " << idTin << " in pltfile to " << idTlocal << " local, to " << idTout << " in output" << endl;
	std::cout << "Y from " << idYin << " in pltfile to " << idYlocal << " local, to " << idCout << " in output" << endl;
    }
    
    PArray<MultiFab> outdata(Nlev,PArrayManage);
    const int nGrow = 0;

    FArrayBox Yt;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        outdata.set(lev,new MultiFab(ba,nCompOut,nGrow));
        MultiFab indata(ba,nCompIn,nGrow);

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        amrData.FillVar(indata,lev,inNames,destFillComps);
        for (int i=0; i<inNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        for (MFIter mfi(indata); mfi.isValid(); ++mfi)
        {
            const FArrayBox& Ylocal = indata[mfi];
            const FArrayBox& Tlocal = indata[mfi];
            FArrayBox& Cout = outdata[lev][mfi];
            FArrayBox& Tout = outdata[lev][mfi];
            const Box& box = mfi.validbox();

            Yt.resize(box, nSpec);
	    int idYtmp = 0;
            cd.normalizeMassFrac(Yt,Ylocal,"N2",box,idYlocal,0);
            cd.massFracToMolarConc (Cout, Yt, Tlocal, Patm, box, idYtmp, idTlocal, idCout);
            Tout.copy(Tlocal,idTlocal,idTout,1);
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Derive finished for level " << lev << endl;
    }

    std::string nfile(getFileRoot(plotFileName) + "_C");

    if (ParallelDescriptor::IOProcessor())
        cout << "Writing new data to " << nfile << endl;
    
    const bool verb = false;
    const AmrData& a = amrData;
    writePlotfile(outdata,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),nfile,outNames,verb);

    BoxLib::Finalize();
    return 0;
}