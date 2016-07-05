#ifndef _COMPCONVECTERMS_F_H_
#define _COMPCONVECTERMS_F_H_ 
#ifndef _COMPDIFFTERMS_F_H_
#define _COMPDIFFTERMS_F_H_ 
#ifndef _COMPCDDRTERMS_F_H_
#define _COMPCDDRTERMS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_PUSHVTOG pushvtog
#    define FORT_COMPCONVECTERMS  compconvecterms
#    define FORT_COMPDIFFTERMS  compdiffterms
#    define FORT_COMPCDDRTERMS  compcddrterms
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_PUSHVTOG       PUSHVTOG
#    define FORT_COMPCONVECTERMS  COMPCONVECTERMS
#    define FORT_COMPDIFFTERMS  COMPDIFFTERMS
#    define FORT_COMPCDDRTERMS  COMPCDDRTERMS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_PUSHVTOG       pushvtog
#    define FORT_COMPCONVECTERMS  compconvecterms
#    define FORT_COMPDIFFTERMS  compdiffterms
#    define FORT_COMPCDDRTERMS  compcddrterms
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_PUSHVTOG       pushvtog_
#    define FORT_COMPCONVECTERMS  compconvecterms_
#    define FORT_COMPDIFFTERMS  compdiffterms_
#    define FORT_COMPCDDRTERMS  compcddrterms_
#  endif

#include <ArrayLim.H>

extern "C" {
   
   void FORT_PUSHVTOG(const int* lo,  const int* hi,
                      const int* dlo, const int* dhi,
                      Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                      const int* nc);
   
   void FORT_COMPCONVECTERMS(const int* lo, const int* hi,
                             const Real* Rho, ARLIM_P(Rholo), ARLIM_P(Rhohi),
                             const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                             const Real* H, ARLIM_P(Hlo), ARLIM_P(Hhi),
                             const Real* Ht, ARLIM_P(Htlo), ARLIM_P(Hthi),
                             const Real* Vx, ARLIM_P(Vxlo), ARLIM_P(Vxhi),
                             const Real* Vy, ARLIM_P(Vylo), ARLIM_P(Vyhi),
                             const Real* C, ARLIM_P(Clo), ARLIM_P(Chi),
                             const Real* dx, const int* ns);
   
   void FORT_COMPDIFFTERMS(const int* lo, const int* hi,
                           const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                           const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                           const Real* H, ARLIM_P(Hlo), ARLIM_P(Hhi),
                           const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
                           const Real* L, ARLIM_P(Llo), ARLIM_P(Lhi),
                           const Real* D, ARLIM_P(Dlo), ARLIM_P(Dhi),
                           const Real* dx, const int* ns);
   
   void FORT_COMPCDDRTERMS(const int* lo, const int* hi,
                           const Real* Rho, ARLIM_P(Rholo), ARLIM_P(Rhohi),
                           const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                           const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                           const Real* Cpmix, ARLIM_P(Cpmixlo), ARLIM_P(Cpmixhi),
                           const Real* H, ARLIM_P(Hlo), ARLIM_P(Hhi),
                           const Real* Ht, ARLIM_P(Htlo), ARLIM_P(Hthi),
                           const Real* Vx, ARLIM_P(Vxlo), ARLIM_P(Vxhi),
                           const Real* Vy, ARLIM_P(Vylo), ARLIM_P(Vyhi),
                           const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
                           const Real* L, ARLIM_P(Llo), ARLIM_P(Lhi),
                           const Real* C, ARLIM_P(Clo), ARLIM_P(Chi),
                           const Real* dx, const int* ns);
   
}


#endif 
#endif
#endif
#endif


