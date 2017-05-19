
#ifndef _GRADT_F_H_
#define _GRADT_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_PUSHVTOG pushvtog
#    define FORT_GRAD     grad
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_PUSHVTOG PUSHVTOG
#    define FORT_GRAD     GRAD
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_PUSHVTOG pushvtog
#    define FORT_GRAD     grad
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_PUSHVTOG pushvtog_
#    define FORT_GRAD     grad_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_PUSHVTOG(const int* lo,  const int* hi,
                       const int* dlo, const int* dhi,
                       Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                       const int* nc);

    void FORT_GRAD(const int* lo, const int* hi,
                   Real*       i, ARLIM_P(ilo), ARLIM_P(ihi),
                   const Real* o, ARLIM_P(olo), ARLIM_P(ohi),
                   const Real* delta);
}
#endif

#endif /*_GRADT_F_H_*/
