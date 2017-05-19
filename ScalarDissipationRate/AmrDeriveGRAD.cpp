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

#ifndef WIN32
#include <unistd.h>
#endif

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "FArrayBox.H"
#include "Geometry.H"
#include "AmrDeriveGRAD_F.H"
#include "AmrDeriveUtilities.H"
#include "ChemDriver.H"
static
void
printErrorExit (std::string message)
{
    if (ParallelDescriptor::IOProcessor ())
    {
	std::cout 
	    << "\nerror: " << message
	    << "\n"
	    << std::endl;
    }
    BoxLib::Finalize ();
    exit(0);
}

static
void
printUsage (char* progName)
{
    if (ParallelDescriptor::IOProcessor ())
    {
	std::cout 
	    << "\n   PROGRAM " << progName
	    << "\n"
	    << "\n        form the gradient of a component in a plotfile.  the component"
            << "\n        and gradient are written to a new plotfile with suffix _grad"
            << "\n        unless append"
	    << "\n"
	    << "\n   required keywords:"
	    << "\n"
	    << "\n       infile = name of plotfile to be read"
 	    << "\n"
	    << "\n   optional keywords:"
	    << "\n"
	    << "\n       append = true or false append derivatives (default false)"
	    << "\n    component = name of component to differentiate (default temp)"
	    << "\n         help = true or false (default false)"
	    << "\n        level = finest level to differentiate (default all)"
	    << "\n       MFbase = base name for appended multifabs (deault NEWDAT)"
	    << "\n      verbose = >0 progress report, >1 AmrData output (default 0)"
	    << "\n"
	    << "\n   required geometry keywords:"
	    << "\n"
	    << "\n       geometry.coord_sys = 0 (Cartesian) or (1 cylindrical)"
	    << "\n       geometry.prob_lo   = low ends of domain"
	    << "\n       geometry.prob_hi   = high ends of domain"
	    << "\n"
	    << "\n   optional geometry keywords:"
	    << "\n"
	    << "\n       geometry.is_periodic = 0 or 1 for each axis (default 0)"
	    << "\n"
	    << std::endl;
    }
    BoxLib::Finalize ();
    exit(0);
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
            
//            if (!BoxLib::UtilCreateDirectory(Level, 0755))
//                BoxLib::CreateDirectoryFailed(Level);
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
main (int argc, char* argv[])
{   
    BoxLib::Initialize (argc, argv);

    //
    // read parameters
    //

    ParmParse pp;

    bool help = (argc < 2);
    pp.query ("help", help);
    if (help)
        printUsage (argv[0]);

    int verbose = 0; 
    pp.query ("verbose", verbose);
    if (verbose > 1)
        AmrData::SetVerbose (true);

    std::string infile; 
    pp.get ("infile", infile);

    std::string component = "Z"; 
    pp.query ("component", component);

    //
    // process the input file
    //
    Real Patm=60; pp.query("Patm",Patm);
    Real P ;
        DataServices::SetBatchMode ();
    Amrvis::FileType fileType (Amrvis::NEWPLT);

    DataServices dataServices(infile, fileType);
    if ( ! dataServices.AmrDataOk()) 
    {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel ();
    pp.query ("finestLevel", finestLevel);
    int number_of_levels = finestLevel + 1;

    bool found = false;
    int component_index_in_file;
    const Array <std::string>& plotVarNames = amrData.PlotVarNames ();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
        if (plotVarNames[i] == component) 
	{
	    component_index_in_file = i;
	    found = true;
	}
    }
    if (!found)
	printErrorExit ("cannot find component " + component + " in plotfile " + infile);

    const int component_index_locally = 0;
    const int number_of_components_to_read = component_index_locally + 1;
    ChemDriver cd;
     
    

    Array <std::string> inVarNames (number_of_components_to_read);
    inVarNames[0] = plotVarNames[component_index_in_file];
    Array <int> destFillComps (number_of_components_to_read);
    destFillComps[0] = 0;

    const int index_of_new_data = number_of_components_to_read;
    const int number_of_components_to_write = index_of_new_data   +1;
    const  int nSpec = cd.numSpecies();
    Array<std::string> speciesN(nSpec+2);
    Array<int> Comp(nSpec+2);
    for (int specCount=0;specCount <nSpec; ++specCount) {
      speciesN[specCount] ="Y("+ cd.speciesNames()[specCount]+")";
      Comp[specCount] =specCount;
   //    if(ParallelDescriptor::IOProcessor())
   //       std::cout << speciesN[specCount]<<std::endl;
                 } 
    speciesN[nSpec] = "temp";
    speciesN[nSpec+1] = "density";
    Comp[nSpec]=nSpec;
    Comp[nSpec+1]=nSpec+1;
    PArray<MultiFab> forRhoD(number_of_levels,PArrayManage);//what is Parraymanage
    
    PArray <MultiFab> local_container (number_of_levels, PArrayManage);
    PArray <Geometry> geoms (number_of_levels, PArrayManage);
    const int nGrow = 2; // for fourth order derivative

    FArrayBox tmp;
    FArrayBox CP;
    for (int lev=0; lev<number_of_levels; ++lev)
    {
        const BoxArray ba = amrData.boxArray(lev);
        forRhoD.set(lev, new MultiFab(ba,2*nSpec+3,0));
        forRhoD[lev].setVal(0);
        amrData.FillVar(forRhoD[lev],lev,speciesN,Comp);
        local_container.set (lev, new MultiFab (ba, number_of_components_to_write, nGrow));
        const Array <Real>& delta = amrData.DxLevel ()[lev];

        // Get input local_container data onto intersection ba
        const int myNComp = 1; // gonna need this for fortran calls

        if (verbose > 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Reading data for level " << lev << std::endl;

        // Initialize grow cells to 0....here don't care to get grads right on c-f
        local_container[lev].setBndry (0.0); 
        amrData.FillVar (local_container[lev], lev, inVarNames, destFillComps);
        for (int i=0; i<inVarNames.size(); ++i)
            amrData.FlushGrids (amrData.StateNumber (inVarNames[i]));

        if (verbose > 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Data has been read for level " << lev << std::endl;

        const bool do_corners = false;
        geoms.set (lev, new Geometry (amrData.ProbDomain ()[lev]));

        // Fix up grow cells.  Use extrap for guess
        const Box& dbox = amrData.ProbDomain ()[lev];

	for (MFIter mfi (local_container[lev]); mfi.isValid (); ++mfi)
        {
            FArrayBox& fab = local_container[lev][mfi];
            const Box& box = mfi.validbox ();
            FORT_PUSHVTOG
		(box.loVect (),
		 box.hiVect (),
		 dbox.loVect (),
		 dbox.hiVect (),
		 fab.dataPtr (component_index_locally),
		 ARLIM (fab.loVect ()),
		 ARLIM (fab.hiVect ()),
		 &myNComp);
		 }

        // Fix up fine-fine and periodic for component_index_locally
        local_container[lev].FillBoundary (component_index_locally, 1);
//        forRhoD[lev].FillBoundary (component_index_locally, 1);

        //geoms[lev].FillPeriodicBoundary (local_container[lev], component_index_locally, 1, do_corners);

        // Compute gradient.  Result in local_container, comp=index_of_new_data
        FArrayBox cp, rhoY;
        MultiFab MF(ba,1,0);
        for (MFIter mfi(local_container[lev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = local_container[lev][mfi];
            const Box& box = mfi.validbox();
            FArrayBox& rhoD = forRhoD[lev][mfi];
            FArrayBox& mf = MF[mfi];
            cp.resize(box,1);
            rhoY.resize(box,nSpec);
            rhoY.copy(rhoD,0,0,nSpec);
            for (int i=0 ; i < nSpec; i++)
            rhoY.mult(rhoD,nSpec+1,i,1);
           cd.getCpmixGivenTY(cp,rhoD,rhoD,box,nSpec,0,0);
           cd.getMixAveragedRhoDiff(rhoD, rhoY, rhoD, Patm, box,0,nSpec,nSpec+2);           

// take thermal diffusivity as the average mass diffusivity

/*           rhoD.copy(rhoD,nSpec+2,2*nSpec+2,1);
           rhoD.mult(rhoD,0,2*nSpec,1);
           for (int is =1 ; is <nSpec ; is++ )
               {  
                  mf.copy(rhoD,nSpec+2+is,0,1);
                  mf.mult(rhoD,is,0,1);
                  rhoD.plus(mf,0,2*nSpec+2,1);
                
                       }
*/
           rhoD.divide(rhoD,nSpec+1,2*nSpec+2,1);   
           rhoD.divide(cp,0,2*nSpec+2,1);
           rhoD.mult(2.0,2*nSpec+2,1);
                  
          
                                   
            FORT_GRAD (box.loVect (),
		       box.hiVect (),
		       fab.dataPtr (component_index_locally),
		       ARLIM (fab.loVect ()),
		       ARLIM (fab.hiVect ()),
		       fab.dataPtr (index_of_new_data), 
		       ARLIM (fab.loVect ()),
		       ARLIM (fab.hiVect ()),
		       delta.dataPtr ());
       
          fab.mult(rhoD,2*nSpec+2,1,1);
          mf.copy(rhoD,2*nSpec+2,0,1); 
          
                       
        }
//       if(lev==1)
//        VisMF::Write(MF,"MF");
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Gradient has been computed on level " << lev << std::endl;
    }

    Array <std::string> names_of_local_components(number_of_components_to_write);
    for (int i=0; i<number_of_components_to_read; ++i)
        names_of_local_components[i] = inVarNames[i];
    names_of_local_components[index_of_new_data+0] = "Scalar_diss";
//    names_of_local_components[index_of_new_data+1] = component + "_gy";
//    names_of_local_components[index_of_new_data+2] = component + "_mag";
//#if BL_SPACEDIM==3
//    names_of_local_components[index_of_new_data+2] = component + "_gz";
//    names_of_local_components[index_of_new_data+3] = component + "_mag";
//#endif

    bool append = true; 
    pp.query ("append", append);
    if (append)
    {
        int number_of_components_to_append = number_of_components_to_write - number_of_components_to_read;
        PArray<MultiFab> container_to_write(number_of_levels, PArrayManage);
        for (int lev=0; lev<number_of_levels; ++lev)
        {
            const BoxArray ba = local_container[lev].boxArray ();
            container_to_write.set (lev, new MultiFab (ba, number_of_components_to_append, nGrow));
            MultiFab::Copy (container_to_write[lev], 
			    local_container[lev], 
			    number_of_components_to_read, 
			    0, 
			    number_of_components_to_append, 
			    0);
        }
        Array <std::string> names_to_write (number_of_components_to_append);
        for (int i=0; i<number_of_components_to_append; ++i)
            names_to_write[i] = names_of_local_components[number_of_components_to_read+i];

        std::string MFbase = "Scalar_diss"; 
	pp.query ("MFbase", MFbase);
        appendToPlotFile (amrData, container_to_write, infile, names_to_write, MFbase,1);
        
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
        {
            std::cout << "...finished.  Note: to see new data, you must rename NewHeader in the" << std::endl;
            std::cout << "              pltfile to Header (probably want to save the original first)" << std::endl;
        }
    }
    else
    {
        std::string outfile (getFileRoot (infile) + "_grad");
        pp.query("outfile",outfile);
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Writing new data to " << outfile << std::endl;

        writePlotfile
	    (local_container,
	     amrData.Time (),
	     amrData.ProbLo (),
	     amrData.ProbHi (),
	     amrData.RefRatio (),
	     amrData.ProbDomain (),
	     amrData.DxLevel (),
	     amrData.CoordSys (),
	     outfile,
	     names_of_local_components);
    }

    //
    // end
    //

    BoxLib::Finalize();
    return 0;
}
