#include <ArrayLim.H>

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
	       const Array<std::string>&  names);
