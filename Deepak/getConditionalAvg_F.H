#ifndef _GETPDF_F_H_
#define _GETPDF_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_GETPDF getpdf
#    define FORT_GETJPDF getjpdf
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_GETPDF GETPDF
#    define FORT_GETJPDF GETJPDF
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_GETPDF getpdf
#    define FORT_GETJPDF getjpdf
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_GETPDF getpdf_
#    define FORT_GETJPDF getjpdf_
#  endif

#include <ArrayLim.H>

extern "C"
{
   void FORT_GETPDF(const int* lo, const int* hi,
                    const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                    const Real* min, const Real* max,
                    const Real* delta, const int* N);
   void FORT_GETJPDF(const int* lo, const int* hi,
                     const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                     const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                     Real* Bins, ARLIM_P(Binslo), ARLIM_P(Binshi),
                     const Real* minX, const Real* maxX,
                     const Real* deltaX, const int* Nx, 
                     const Real* minY, const Real* maxY,
                     const Real* deltaY, const int* Ny, const Real* vol,
                     const Real* mask, ARLIM_P(masklo),ARLIM_P(maskhi));
} 

#endif
#endif
