#ifndef _GETSUMRR_F_H_
#define _GETSUMRR_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_GETSUMRR getsumrr
#    define FORT_GETSUMHEATRELEASE getsumheatrelease
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_GETSUMRR GETSUMRR
#    define FORT_GETSUMHEATRELEASE GETSUMHEATRELEASE
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_GETSUMRR getsumrr
#    define FORT_GETSUMHEATRELEASE getsumheatrelease
#  elif defined(BL_FORT_USE_UNDERSCORE)
#     define FORT_GETSUMRR getsumrr_
#    define FORT_GETSUMHEATRELEASE getsumheatrelease_
#  endif


#include <REAL.H>
#include <ArrayLim.H>

extern "C"
{
   void FORT_GETSUMRR(const int* lo, const int* hi,
                      const Real* FWD, ARLIM_P(FWDlo), ARLIM_P(FWDhi),
                      const Real* REV, ARLIM_P(REVlo), ARLIM_P(REVhi),
                      const Real* SUM, ARLIM_P(SUMlo), ARLIM_P(SUMhi),
                      const int* Nreacs);
   void FORT_GETSUMHEATRELEASE(const int* lo, const int* hi,
                      const Real* HR, ARLIM_P(HRlo), ARLIM_P(HRhi),
                                const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi), 
                                const Real* minT, const Real* maxT,
                               const Real* deltaT, const int* Ns, const Real* Sum);
} 

#endif
#endif
