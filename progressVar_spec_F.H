
#ifndef _COMPDIFFTERMS_F_H_
#define _COMPDIFFTERMS_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_PUSHVTOG pushvtog
#    define FORT_COMPDIFFTERMS  compdiffterms
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_PUSHVTOG       PUSHVTOG
#    define FORT_COMPDIFFTERMS  COMPDIFFTERMS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_PUSHVTOG       pushvtog
#    define FORT_COMPDIFFTERMS  compdiffterms
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_PUSHVTOG       pushvtog_
#    define FORT_COMPDIFFTERMS  compdiffterms_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_PUSHVTOG(const int* lo,  const int* hi,
                       const int* dlo, const int* dhi,
                       Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                       const int* nc);

    void FORT_COMPDIFFTERMS(const int* lo, const int* hi,
                            const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                            Real* out,
                             Real* Max);
}
#endif

#endif /*_COMPDIFFTERMS_F_H_*/
