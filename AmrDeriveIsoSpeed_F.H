#ifndef _AMRDERIVEISO_F_H_
#define _AMRDERIVEISO_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_SETLOC  setloc
#    define FORT_SETCLOC  setcloc
#    define FORT_COMPDIFFTERMS compdiffterms
#    define FORT_PUSHVTOG pushvtog
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_SETLOC  SETLOC
#    define FORT_SETCLOC SETCLOC
#    define FORT_PUSHVTOG PUSHVTOG
#    define FORT_COMPDIFFTERMS COMPDIFFTERMS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_SETLOC  setloc
#    define FORT_SETCLOC setcloc
#    define FORT_COMPDIFFTERMS compdiffterms
#    define FORT_PUSHVTOG pushvtog
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_TRICUBIC_GET_COEFF_STACKED fort_tricubic_get_coeff_stacked_
#    define FORT_TRICUBIC_GET_COEFF fort_tricubic_get_coeff_
#    define FORT_TRICUBIC_EVAL fort_tricubic_eval_
#    define FORT_SETLOC  setloc_   
#    define FORT_SETCLOC setcloc_
#    define FORT_COMPDIFFTERMS compdiffterms_
#    define FORT_PUSHVTOG pushvtog_
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
                            const Real* dx, int* indices, std::vector<Array<Real>, std::allocator<Array<Real> > >* GridPts,  const int* nNodes,
                            const Real* data, ARLIM_P(dataLo),ARLIM_P(dataHi),const int* auxcomps, const int* tot_comps, const int* sz);
  void FORT_PUSHVTOG(const int* lo,  const int* hi,
                       const int* dlo, const int* dhi,
                       Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                       const int* nc);


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
