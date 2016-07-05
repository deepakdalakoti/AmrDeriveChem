#ifndef _CONDITIONEDRR_F_H_
#define _CONDITIONEDRR_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_CONDITIONEDRR CONDITIONEDRR
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CONDITIONEDRR CONDITIONEDRR
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CONDITIONEDRR conditionedrr
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CONDITIONEDRR conditionedrr_
#  endif

#include <ArrayLim.H>

extern "C"
{
   void FORT_CONDITIONEDRR(const int* lo, const int* hi,
                           const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                           const Real* Fw, ARLIM_P(Fwlo), ARLIM_P(Fwhi),
                           const Real* Bw, ARLIM_P(Bwlo), ARLIM_P(Bwhi),
                           const Real* lowval, const Real* upval, const int*Nreacs);
} 

#endif
#endif
