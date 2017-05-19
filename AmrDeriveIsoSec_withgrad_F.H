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
#    define FORT_TRICUBIC_GET_COEFF_STACKED fort_tricubic_get_coeff_stacked_
#    define FORT_TRICUBIC_GET_COEFF fort_tricubic_get_coeff_
#    define FORT_TRICUBIC_EVAL fort_tricubic_eval_
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
    void FORT_COMPDIFFTERMS(const int* lo, const int* hi,const int* flo, const int* fhi,
                            const int* slo,const int* shi,ARLIM_P(flo), ARLIM_P(shi),
                            ARLIM_P(slo), ARLIM_P(shi),
                            const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                            const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
                            const Real* Reac, ARLIM_P(Reaclo),
                            ARLIM_P(Reachi),const Real* density, ARLIM_P(densitylo), 
                            ARLIM_P(densityhi),const Real* T, Real* Speed,const Real* dx, const int* ws, const int* ns,  const Real* velocity,
                            const Real* pos, const Real* Z, const Real* alpha);


   void FORT_TRICUBIC_GET_COEFF_STACKED(Real*       a,
                                         const Real* x);

    void FORT_TRICUBIC_GET_COEFF(Real*       a,
                                 const Real* f,
                                 const Real* dfdx,
                                 const Real* dfdy,
                                 const Real* dfdz,
                                 const Real* d2fdxdy,
                                 const Real* d2fdxdz,
                                 const Real* d2fdydz,
                                 const Real* d3fdxdydz);
    
    void FORT_TRICUBIC_EVAL(const Real* a,
                            const Real* x,
                            const Real* y,
                            const Real* z,
                            Real*       eval);



}
#endif

#endif /*_AMRDERIVEJUNK_F_H_*/
