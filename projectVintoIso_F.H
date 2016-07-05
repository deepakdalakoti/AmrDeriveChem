
#ifndef _PROJECTVINTOISO_F_H_
#define _PROJECTVINTOISO_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_PUSHVTOG pushvtog
#    define FORT_PRJECTV  prjectv
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_PUSHVTOG PUSHVTOG
#    define FORT_PRJECTV  PRJECTV
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_PUSHVTOG pushvtog
#    define FORT_PRJECTV  prjectv
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_PUSHVTOG pushvtog_
#    define FORT_PRJECTV  prjectv_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_PUSHVTOG(const int* lo,  const int* hi,
                       const int* dlo, const int* dhi,
                       Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                       const int* nc);

    void FORT_PRJECTV(const int* lo, const int* hi,
                      const Real* i,  ARLIM_P(ilo),  ARLIM_P(ihi),
                      const Real* ii, ARLIM_P(iilo), ARLIM_P(iihi),
                      const Real* o,  ARLIM_P(olo),  ARLIM_P(ohi),
                      const Real* delta);
}
#endif

#endif /*_PROJECTVINTOISO_F_H_*/
