
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

    void FORT_COMPDIFFTERMS(const int* flo, const int* fhi,const int* lo, const int* hi,
                            const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                            const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
                            const Real* density, ARLIM_P(densitylo),
                            ARLIM_P(densityhi), const Real* vel, ARLIM_P(vello),
                            ARLIM_P(velhi),const Real* D, ARLIM_P(Dlo),
                            ARLIM_P(Dhi),const Real* dx, const int* ns, const int* si,
                            const int* nwhich, Real* Adv, Real* total,const Real* Ncomp);
}
#endif

#endif /*_COMPDIFFTERMS_F_H_*/
