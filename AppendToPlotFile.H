// --------------------------------------------------------------------
// AppendToPlotFile.H
// --------------------------------------------------------------------
#ifndef BL_APPENDTOPLOTFILE_H
#define BL_APPENDTOPLOTFILE_H

#include <REAL.H>
#include <FArrayBox.H>
#include <ParallelDescriptor.H>
#include <DataServices.H>
#include <Utility.H>
#include <PArray.H>

void
AppendToPlotFile(AmrData&                  amrData,
                 const PArray<MultiFab>&   mfout,
                 std::string&              oFile,
                 const Array<std::string>& nnames,
                 const std::string&        mfBaseName,
                 const std::string&        NewHeaderName,
                 bool                      verbose);

#endif
// --------------------------------------------------------------------
