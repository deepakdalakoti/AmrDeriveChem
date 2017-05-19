// --------------------------------------------------------------------
// WritePlotFile.H
// --------------------------------------------------------------------
#ifndef BL_WRITEPLOTFILE_H
#define BL_WRITEPLOTFILE_H

#include <REAL.H>
#include <FArrayBox.H>
#include <ParallelDescriptor.H>
#include <DataServices.H>
#include <PArray.H>
#include <MultiFab.H>
#include <Utility.H>


void WritePlotfile(const std::string         &pfversion,
                   const PArray<MultiFab>    &data,
                   const Real                 time,
                   const Array<Real>         &probLo,
                   const Array<Real>         &probHi,
                   const Array<int>          &refRatio,
                   const Array<Box>          &probDomain,
                   const Array<Array<Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const Array<std::string>  &names,
                   const bool                 verbose,
		   const bool                 isCartGrid = false,
		   const Real                *vfeps = 0,
		   const int                 *levelSteps = 0);



// this function makes a 2d plotfile from a 3d one
// the 3d file must have width = 1 in one direction
// compile in 3d

void WritePlotfile2DFrom3D(const std::string         &pfversion,
                           const PArray<MultiFab>    &data,
                           const Real                 time,
                           const Array<Real>         &probLo,
                           const Array<Real>         &probHi,
                           const Array<int>          &refRatio,
                           const Array<Box>          &probDomain,
                           const Array<Array<Real> > &dxLevel,
                           const int                  coordSys,
                           const std::string         &oFile,
                           const Array<std::string>  &names,
                           const bool                 verbose,
		           const bool                 isCartGrid = false,
		           const Real                *vfeps = 0,
		           const int                 *levelSteps = 0);

#endif
// --------------------------------------------------------------------
