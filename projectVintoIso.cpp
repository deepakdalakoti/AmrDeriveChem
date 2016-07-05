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
#include "projectVintoIso_F.H"

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
    
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if (verbose && ParallelDescriptor::IOProcessor())
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
    const int finestLevel = data.size() - 1;
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
    
    if (verbose && ParallelDescriptor::IOProcessor())
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

    int finestLevel = mfout.size() - 1;
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

    bool verbose=false;
    if (pp.contains("verbose"))
    {
        verbose=true;
        AmrData::SetVerbose(true);
    }

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    std::string progressName = "temp"; pp.query("progressName",progressName);
    Array<int> vComps(BL_SPACEDIM,-1);
    int nc = pp.countval("vel");
    if (nc!=BL_SPACEDIM)
        BoxLib::Abort("must specify a velocity component for each spatial dimension!");
    pp.getarr("vel",vComps,0,nc);

    int idC = -1;
    const Array<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == progressName) idC = i;
    }
    if (ParallelDescriptor::IOProcessor())
    {
        if (idC<0)
            cerr << "Cannot find required data in pltfile" << endl;
        for (int i=0; i<BL_SPACEDIM; ++i)
            if (vComps[i]<0)
                cerr << "Cannot find vector component in pltfile" << endl;
    }

    const int idCst = 0;
    const int idVst = idCst + 1;
    const int nCompIn = idVst + BL_SPACEDIM;

    Array<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];
    for (int i=0; i<BL_SPACEDIM; ++i)
        inVarNames[idVst+i] = plotVarNames[vComps[i]];
    if (ParallelDescriptor::IOProcessor() && verbose)
        for (int i=0; i<inVarNames.size(); ++i)
            std::cout << "going to read: " << inVarNames[i] << endl;

    Array<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i)
        destFillComps[i] = i;

    const int idPV = nCompIn;
    const int nCompOut = idPV + BL_SPACEDIM;

    PArray<MultiFab> state(Nlev,PArrayManage);
    PArray<Geometry> geoms(Nlev,PArrayManage);
    const int nGrow = 1;

    FArrayBox tmp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        state.set(lev,new MultiFab(ba,nCompOut,nGrow));
        const Array<Real>& delta = amrData.DxLevel()[lev];

        if (ParallelDescriptor::IOProcessor())
            cerr << "Reading data for level " << lev << endl;

        state[lev].setBndry(0.0); // Initialize grow cells to 0....here don't care to get grads right on c-f
        amrData.FillVar(state[lev],lev,inVarNames,destFillComps);
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));

        if (ParallelDescriptor::IOProcessor())
            cerr << "Data has been read for level " << lev << endl;

        const bool do_corners = false;
        geoms.set(lev,new Geometry(amrData.ProbDomain()[lev]));

        // Fix up grow cells.  Use extrap for guess
        const Box& dbox = amrData.ProbDomain()[lev];
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_PUSHVTOG(box.loVect(),box.hiVect(),dbox.loVect(),dbox.hiVect(),
                          s.dataPtr(idCst),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                          &nCompIn);
        }

        // Fix up fine-fine and periodic for idCst
        state[lev].FillBoundary(idCst,nCompIn);
        geoms[lev].FillPeriodicBoundary(state[lev],idCst,nCompIn,do_corners);

        // Project vector into progress contours
        for (MFIter mfi(state[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = state[lev][mfi];
            const Box& box = mfi.validbox();
            FORT_PRJECTV(box.loVect(),box.hiVect(),
                         s.dataPtr(idCst), ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         s.dataPtr(idVst), ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         s.dataPtr(idPV),  ARLIM(s.loVect()),ARLIM(s.hiVect()),
                         delta.dataPtr());
        }

        if (ParallelDescriptor::IOProcessor())
            cerr << "Vector has been projected on level " << lev << endl;

    }

    Array<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i)
        nnames[i] = inVarNames[i];
    nnames[idPV+0] = plotVarNames[vComps[0]] + "_pgx";
    nnames[idPV+1] = plotVarNames[vComps[1]] + "_pgy";
#if BL_SPACEDIM==3
    nnames[idPV+2] = plotVarNames[vComps[2]] + "_pgz";
#endif

    bool verb=false;
    bool appendPlotFile=false; pp.query("appendPlotFile",appendPlotFile);
    if (appendPlotFile)
    {

        int nStateOut = nCompOut - nCompIn;
        PArray<MultiFab> ostate(Nlev,PArrayManage);
        for (int lev=0; lev<Nlev; ++lev)
        {
            const BoxArray ba = state[lev].boxArray();
            ostate.set(lev,new MultiFab(ba,nStateOut,nGrow));
            MultiFab::Copy(ostate[lev],state[lev],nCompIn,0,nStateOut,0);
        }
        Array<std::string> namesOut(nStateOut);
        for (int i=0; i<nStateOut; ++i)
            namesOut[i] = nnames[nCompIn+i];

        std::string newMFBaseName = "NEWDAT"; pp.query("newMFBaseName",newMFBaseName);
        appendToPlotFile(amrData,ostate,plotFileName,namesOut,newMFBaseName,verb);
        
        if (ParallelDescriptor::IOProcessor())
        {
            cerr << "...finished.  Note: to see new data, you must rename NewHeader in the" << endl;
            cerr << "              pltfile to Header (probably want to save the original first)" << endl;
        }
    }
    else
    {
        std::string nfile(getFileRoot(plotFileName) + "_proj");

        if (ParallelDescriptor::IOProcessor())
            cout << "Writing new data to " << nfile << endl;

        const AmrData& a = amrData;
        writePlotfile(state,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                      a.DxLevel(),a.CoordSys(),nfile,nnames,verb);
    }

    BoxLib::Finalize();
    return 0;
}
