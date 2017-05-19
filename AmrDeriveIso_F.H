#ifndef _AMRDERIVEISO_F_H_
#define _AMRDERIVEISO_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_SETLOC  setloc
#    define FORT_SETCLOC  setcloc
#    define FORT_COMPDIFFTERMS compdiffterms
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_SETLOC  SETLOC
#    define FORT_SETCLOC SETCLOC
#    define FORT_COMPDIFFTERMS COMPDIFFTERMS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_SETLOC  setloc
#    define FORT_SETCLOC setcloc
#    define FORT_COMPDIFFTERMS compdiffterms
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_SETLOC  setloc_
#    define FORT_SETCLOC setcloc_
#    define FORT_COMPDIFFTERMS compdiffterms_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_SETLOC(const int* lo,  const int* hi,
                     Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                     const Real* dx, const Real* domlo);

    void FORT_SETCLOC(const int* lo,  const int* hi,
                     Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                     const Real* fdx, const Real* domnlo, const int* ratio);
    void FORT_COMPDIFFTERMS(const int* lo, const int* hi,
                            const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                            const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
                            const Real* Reac, ARLIM_P(Reaclo),
                            ARLIM_P(Reachi),const Real* density, ARLIM_P(densitylo),
                            ARLIM_P(densityhi), Real* Speed,const Real* dx, const int* ws, const int* ns,  const Real* velocity,
                            const Real* pos, const Real* Z);
}
#endif

#endif /*_AMRDERIVEJUNK_F_H_*/
