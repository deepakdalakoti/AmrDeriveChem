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

static
vector<std::string>
Tokenize (const std::string& instr,
          const std::string& separators)
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

void WritePlotFile(const PArray<MultiFab>&   mfa,
                   const Array<std::string>& names,
		   AmrData&                  amrdToMimic,
		   const std::string&        oFile,
                   const Box&                subbox,
		   bool                      verbose);

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
    std::string infile;

    bool verbose = verbose_DEF;
    verbose = (pp.contains("verbose") ? true : false);
    if (verbose)
        AmrData::SetVerbose(true);

    pp.get("infile",infile);

    vector<std::string> pieces = Tokenize(infile,std::string("/"));
    std::string outfile = pieces[pieces.size()-1] + std::string("_section");
    pp.query("outfile",outfile);

    // Read in pltfile and get amrData ref
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = dataServices.AmrDataRef();

    Array<int> comps;
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
        BL_ASSERT(sComp+nComp < amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

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
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        // Build the BoxArray for the result
        BoxList bl;
        const BoxArray& ba_all = amrData.boxArray(iLevel);
        for (int i=0; i<ba_all.size(); ++i)
        {
            const Box ovlp = subboxes[iLevel] & ba_all[i];
            if (ovlp.ok())
                bl.push_back(ovlp);
        }
        BoxArray ba_sub(bl);
        if (ba_sub.size() > 0)
        {
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
    WritePlotFile(data_sub,names,amrData,outfile,subboxes[0],verbose);
    
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

void WritePlotFile(const PArray<MultiFab>& mfa,
                   const Array<std::string>&   names,
		   AmrData&                amrdToMimic,
		   const std::string&          oFile,
                   const Box&              subbox,
		   bool                    verbose)
{
    int ntype = names.size();
    int finestLevel = mfa.size() - 1;
    
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
  
    os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "Opening file = " << oFileHeader << '\n';

#ifdef BL_USE_NEW_HFILES
    os.open(oFileHeader.c_str(), ios::out|ios::binary);
#else
    os.open(oFileHeader.c_str(), ios::out);
#endif

    if (os.fail())
        BoxLib::FileOpenFailed(oFileHeader);
    //
    // Start writing plotfile.
    //
    os << amrdToMimic.PlotFileVersion() << '\n';
    int n_var = ntype;
    os << n_var << '\n';
    for (int n = 0; n < ntype; n++) os << names[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << amrdToMimic.Time() << '\n';
    os << finestLevel << '\n';
    int i;
    const int nLev = finestLevel + 1;
    Real ProbLo[BL_SPACEDIM], ProbHi[BL_SPACEDIM];
    for (i = 0; i < BL_SPACEDIM; ++i)
    {
        ProbLo[i] = amrdToMimic.ProbLo()[i]
            + amrdToMimic.DxLevel()[0][i]
            * (subbox.loVect()[i] - amrdToMimic.ProbDomain()[0].loVect()[i]);
        ProbHi[i] = amrdToMimic.ProbHi()[i]
            - amrdToMimic.DxLevel()[0][i]
            * (amrdToMimic.ProbDomain()[0].hiVect()[i] - subbox.hiVect()[i]);
    }
    Array<Box> ProbDomain(nLev);
    for (i = 0; i < nLev; ++i)
    {
        if (i==0) 
        {
            ProbDomain[i] = subbox;
        }
        else
        {
            ProbDomain[i] = BoxLib::refine(ProbDomain[i-1],amrdToMimic.RefRatio()[i-1]);
        }
    }

    for (i = 0; i < BL_SPACEDIM; i++) os << ProbLo[i] << ' ';
    os << '\n';
    for (i = 0; i < BL_SPACEDIM; i++) os << ProbHi[i] << ' ';
    os << '\n';
    for (i = 0; i < finestLevel; i++) os << amrdToMimic.RefRatio()[i] << ' ';
    os << '\n';
    for (i = 0; i <= finestLevel; i++) os << ProbDomain[i] << ' ';
    os << '\n';
    for (i = 0; i <= finestLevel; i++) os << 0 << ' ';
    os << '\n';
    for (i = 0; i <= finestLevel; i++)
    {
        for (int k = 0; k < BL_SPACEDIM; k++)
            os << amrdToMimic.DxLevel()[i][k] << ' ';
        os << '\n';
    }
    os << amrdToMimic.CoordSys() << '\n';
    os << "0\n"; // The bndry data width.
    //
    // Write out level by level.
    //
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Write state data.
        //
        BoxArray ba_sub = mfa[iLevel].boxArray();
        int nGrids = ba_sub.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
    
        if (ParallelDescriptor::IOProcessor())
        {
            os << iLevel << ' ' << nGrids << ' ' << amrdToMimic.Time() << '\n';
            os << 0 << '\n';
    
            for (i = 0; i < nGrids; ++i)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    os << amrdToMimic.DxLevel()[iLevel][n]*ba_sub[i].smallEnd(n)
                       << ' '
                       << amrdToMimic.DxLevel()[iLevel][n]*(ba_sub[i].bigEnd(n)+1)
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
        VisMF::Write(mfa[iLevel], PathName, VisMF::OneFilePerCPU);
    }

    os.close();
}
