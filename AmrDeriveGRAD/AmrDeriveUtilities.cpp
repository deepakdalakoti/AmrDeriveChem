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

void
writePlotfile (const PArray<MultiFab>&    data,
	       Real                       time,
	       const Array<Real>&         probLo,
	       const Array<Real>&         probHi,
	       const Array<int>&          refRatio,
	       const Array<Box>&          probDomain,
	       const Array<Array<Real> >& dxLevel,
	       int                        coordSys,
	       std::string&               outplot,
	       const Array<std::string>&  names)
{
    const int finestLevel = data.size() - 1;
    static const std::string MultiFabBaseName("/Cell");

    //
    // create the directories and header
    //
    if (ParallelDescriptor::IOProcessor())
    {
	// create the plotfile directory
        if (!BoxLib::UtilCreateDirectory(outplot,0755))
            BoxLib::CreateDirectoryFailed(outplot);
    
	// create the plotfile header
	std::ofstream os;

	std::string oFileHeader(outplot);
	oFileHeader += "/Header";
	
	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
	
	os.open (oFileHeader.c_str (), std::ios::out|std::ios::binary);
	if (os.fail ())
	    BoxLib::FileOpenFailed(oFileHeader);
	
	// write the plotfile header
	std::string plotFileVersion = "NavierStokes-V1.1";
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

	// loop over levels only for io processor at this time
	for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	{
	    // write the level information into the header
	    const BoxArray& ba = data[iLevel].boxArray();
	    int nGrids = ba.size();
	    char buf[64];
	    sprintf(buf, "Level_%d", iLevel);
	    
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

	    // build the full relative pathname of the MultiFab.
            std::string RelativePathName(buf);
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';

	    // build the directory to hold the MultiFabs at this level
	    std::string Level(outplot);
	    Level += '/';
	    Level += buf;
	    if (!BoxLib::UtilCreateDirectory(Level, 0755))
		BoxLib::CreateDirectoryFailed(Level);
	}
	// close the stream to write the header
	os.close();
    }
	
    //
    // force all processors to wait until the directories are built
    //
    ParallelDescriptor::Barrier();

    //
    // Write out level by level
    //
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        std::string PathName(outplot);
        PathName += '/';
        PathName += buf;
        PathName += MultiFabBaseName;
        
        VisMF::Write(data[iLevel], PathName);
    }
    
}

