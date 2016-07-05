
#ifndef _CONVERTDRM19TOH2_F_H_
#define _CONVERTDRM19TOH2_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_CONVERT_DRM19_TO_H2 drm19toh2
#    define FORT_TfromHYpt_new       tfrhyptn
#    define CKHBMS_new               ckhbmsn
#    define CKCPBS_new               ckcpbsn
#    define CKINDX_new               ckindxn
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CONVERT_DRM19_TO_H2 DRM19TOH2
#    define FORT_TfromHYpt_new       TFRHYPTN
#    define CKHBMS_new               CKHBMSN
#    define CKCPBS_new               CKCPBSN
#    define CKINDX_new               CKINDXN
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CONVERT_DRM19_TO_H2 drm19toh2
#    define FORT_TfromHYpt_new       tfrhyptn
#    define CKHBMS_new               ckhbmsn
#    define CKCPBS_new               ckcpbsn
#    define CKINDX_new               ckindxn
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CONVERT_DRM19_TO_H2 drm19toh2_
#    define FORT_TfromHYpt_new       tfrhyptn_
#    define CKHBMS_new               ckhbmsn_
#    define CKCPBS_new               ckcpbsn_
#    define CKINDX_new               ckindxn_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_CONVERT_DRM19_TO_H2(const Real* T, Real* Y);

    int FORT_TfromHYpt_new(Real* T, const Real* Hin, const Real* Yin,
                           const Real* errMAX, const int* NiterMAX, const Real* res);

    void CKHBMS_new(double * T, double * y, int * iwrk, double * rwrk, double * hbms);
    void CKCPBS_new(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
    void CKINDX_new(int * iwrk, double * rwrk, int * e, int * s, int * r, int * f);
}

#endif

#endif /*_CONVERTDRM19TOH2_F_H_*/
