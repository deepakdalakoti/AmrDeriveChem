
#ifndef _HACK_F_H_
#define _HACK_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_HACKVAL  hackval
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_HACKVAL  HACKVAL
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_HACKVAL  hackval
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_HACKVAL  hackval_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_HACKVAL(const int* lo, const int* uhi,
                      Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                      const Real* dx, const Real* plo);

}
#endif

#endif /*_HACK_F_H_*/
