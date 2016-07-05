#ifndef _AMRDERIVEISO_F_H_
#define _AMRDERIVEISO_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_SETLOC  setloc
#    define FORT_SETCLOC  setcloc
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_SETLOC  SETLOC
#    define FORT_SETCLOC SETCLOC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_SETLOC  setloc
#    define FORT_SETCLOC setcloc
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_SETLOC  setloc_
#    define FORT_SETCLOC setcloc_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_SETLOC(const int* lo,  const int* hi,
                     Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                     const Real* dx, const Real* domlo);

    void FORT_SETCLOC(const int* lo,  const int* hi,
                     Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                     const Real* fdx, const Real* domnlo, const int* ratio);
}
#endif

#endif /*_AMRDERIVEJUNK_F_H_*/
