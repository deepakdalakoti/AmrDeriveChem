#ifndef _SUMMASSCYL_F_H_
#define _SUMMASSCYL_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_SUM_IN_RZ  fort_sum_in_rz
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_SUM_IN_RZ FORT_SUM_IN_RZ
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_SUM_IN_RZ fort_sum_in_rz
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_SUM_IN_RZ fort_sum_in_rz_
#  endif

#include <ArrayLim.H>

extern "C" {

    Real FORT_SUM_IN_RZ(const int* lo,  const int* hi,
                        const Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                        const Real* dx, const Real* dlo,
                        const Real* rMax, const Real* zMax);
}
#endif

#endif
