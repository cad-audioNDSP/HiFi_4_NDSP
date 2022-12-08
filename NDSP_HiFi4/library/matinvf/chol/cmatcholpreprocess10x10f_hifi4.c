/* ------------------------------------------------------------------------ */
/* Copyright (c) 2019 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/* IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                     */
/*                                                                          */
/* DSP Library                                                              */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2014-2019 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
#include "NatureDSP_types.h"
#include "common.h"
#include "NatureDSP_Signal_matinv.h"
#include <math.h>
#include "common_fpu.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*-------------------------------------------------------------------------
Preprocessing for Least Square Solutions
The result is matrix Z[NxN], such that
Z = A'*A + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]          matrix A. For the fixed point, the representation is Q(qA)
sigma2          regularization term. For fixed point, the representation 
                should be Q(2*qA-30)
qRA             qR-qA; difference between fixed point representations of 
                decomposition matrix R and original matrix A (for the fixed 
                point API only). Should be equal or less than 0 (typically 
                -2).
Output:
Z[N*N]          matrix Z. For the fixed point, the representation is Q(2*qR-4)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
/* single-matrix API */
void  cmatcholpreprocess10x10f(void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t  sigma2)
{
    int m, n, p;
    xtfloatx2 * restrict pRv0;
    //xtfloatx2 * restrict pRv1;
    xtfloatx2 * restrict pRh0;
    //xtfloatx2 * restrict pRh1;
    const xtfloatx2 * restrict pApm;
    const xtfloatx2 * restrict pApn;
    xtfloatx2 Apn0, Apn1;
    xtfloatx2 Apm0, Apm1;
    xtfloatx2 Acc00, Acc01, Acc10, Acc11, Acc00t, Acc01t, Acc10t, Acc11t;
    xtfloatx2 sigma = sigma2;
    xtfloatx2 sigma_w_zeros;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    sigma_w_zeros = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(AE_SLAI64(AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(sigma)), 32)));

    for (m = 0; m < 10; m += 2)
    {
        Acc00 = Acc11 = sigma_w_zeros;
        Acc01 = Acc10 = Acc00t = Acc01t = Acc10t = Acc11t = XT_CONST_S(0);
        pRv0 = pRh0 = (xtfloatx2*)(R + m * (10 + 1));
        //pRv1 = (xtfloatx2*)(R + m * (4 + 1) + 4);
        //pRh1 = (xtfloatx2*)(R + m * (4 + 1) + 1);
        for (n = m; n < 10; n += 2)
        {
            pApm = (xtfloatx2*)(A + m);
            pApn = (xtfloatx2*)(A + n);
            for (p = 0; p < 10; p++)
            {
                XT_LSX2IP(Apm0, pApm, SZ_CF32);
                XT_LSX2XP(Apm1, pApm, (10 - 1)*SZ_CF32);
                XT_LSX2IP(Apn0, pApn, SZ_CF32);
                XT_LSX2XP(Apn1, pApn, (10 - 1)*SZ_CF32);

                XT_MADDMUX_S(Acc00, Apm0, Apn0, 0);
                XT_MADDMUX_S(Acc00t, Apm0, Apn0, 3);
                XT_MADDMUX_S(Acc01, Apm0, Apn1, 0);
                XT_MADDMUX_S(Acc01t, Apm0, Apn1, 3);
                XT_MADDMUX_S(Acc10, Apm1, Apn0, 0);
                XT_MADDMUX_S(Acc10t, Apm1, Apn0, 3);
                XT_MADDMUX_S(Acc11, Apm1, Apn1, 0);
                XT_MADDMUX_S(Acc11t, Apm1, Apn1, 3);
            }
            Acc00 = XT_ADD_SX2(Acc00, Acc00t);
            Acc01 = XT_ADD_SX2(Acc01, Acc01t);
            Acc10 = XT_ADD_SX2(Acc10, Acc10t);
            Acc11 = XT_ADD_SX2(Acc11, Acc11t);

            XT_SSX2IP(Acc00, pRv0, SZ_CF32);
            XT_SSX2XP(Acc01, pRv0, (10 - 1)*SZ_CF32);
            XT_SSX2IP(Acc10, pRv0, SZ_CF32);
            XT_SSX2XP(Acc11, pRv0, -(10 - 1)*SZ_CF32);
            XT_SSX2IP(XT_CONJC_S(Acc00), pRh0, SZ_CF32);
            XT_SSX2XP(XT_CONJC_S(Acc10), pRh0, (10 - 1)*SZ_CF32);
            XT_SSX2IP(XT_CONJC_S(Acc01), pRh0, SZ_CF32);
            XT_SSX2XP(XT_CONJC_S(Acc11), pRh0, (10 - 1)*SZ_CF32);

            Acc00 = Acc01 = Acc10 = Acc11 = Acc00t = Acc01t = Acc10t = Acc11t = XT_CONST_S(0);
        }
    }
}

size_t  cmatcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#elif(HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*-------------------------------------------------------------------------
Preprocessing for Least Square Solutions
The result is matrix Z[NxN], such that
Z = A'*A + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]          matrix A. For the fixed poiont, the representation is Q(qA)
sigma2          regularization term. For fixed point, the representation
should be Q(2*qA-30)
qRA             qR-qA; difference between fixed point representations of
decomposition matrix R and original matrix A (for the fixed
point API only). Should be equal or less than 0 (typically
-2).
Output:
Z[N*N]          matrix Z. For the fixed poiont, the representation is Q(2*qR-4)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholpreprocess10x10f(void* pScr,
    complex_float *R,
    const complex_float * A,
    const float32_t  sigma2)
{
#if 0
    const float32_t* restrict Af = (const float32_t*)A;
    float32_t* restrict Rf = (float32_t*)R;
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    float32_t re_Apm;
    float32_t im_Apm;
    float32_t re_Apn;
    float32_t im_Apn;
    float32_t re_acc;
    float32_t im_acc;
    for (m = 0; m < 10; m++)
    {
        re_acc = sigma2;
        for (p = 0; p < 10; p++)
        {
            re_Apm = Af[(p * 10 + m) * 2 + 0];
            im_Apm = Af[(p * 10 + m) * 2 + 1];

            re_acc += re_Apm*re_Apm + im_Apm*im_Apm;
        }
        Rf[(m * 10 + m) * 2 + 0] = re_acc;
        Rf[(m * 10 + m) * 2 + 1] = 0.f;

        for (n = m + 1; n < 10; n++)
        {
            re_acc = 0.f;//m == n ? sigma2 : 0.f;
            im_acc = 0.f;

            for (p = 0; p < 10; p++)
            {
                re_Apm = Af[(p * 10 + m) * 2 + 0];
                im_Apm = Af[(p * 10 + m) * 2 + 1];
                re_Apn = Af[(p * 10 + n) * 2 + 0];
                im_Apn = Af[(p * 10 + n) * 2 + 1];

                re_acc += re_Apm*re_Apn + im_Apm*im_Apn;
                im_acc += re_Apm*im_Apn - im_Apm*re_Apn;
            }
            Rf[(m * 10 + n) * 2 + 0] = re_acc;
            Rf[(m * 10 + n) * 2 + 1] = im_acc;
            Rf[(m + n * 10) * 2 + 0] = re_acc;
            Rf[(m + n * 10) * 2 + 1] = -im_acc;
        }
    }
#elif 1
    const float32_t* restrict Af = (const float32_t*)A;
    float32_t* restrict Rf = (float32_t*)R;
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    float32_t re_Apm0, re_Apm1;
    float32_t im_Apm0, im_Apm1;
    float32_t re_Apn0, re_Apn1;
    float32_t im_Apn0, im_Apn1;
    float32_t re_acc00, re_acc01, re_acc10, re_acc11;
    float32_t im_acc00, im_acc01, im_acc10, im_acc11;
    /* processing by squares 2x2 */
    for (m = 0; m < 10; m += 2)
    {
        /* diagonal squares */
        re_acc00 = re_acc11 = sigma2;
        re_acc01 = re_acc10 = 0.f;
        im_acc00 = im_acc01 = im_acc10 = im_acc11 = 0.f;

        for (p = 0; p < 10; p++)
        {
            re_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 0];
            re_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 0];
            im_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 1];
            im_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 1];

            re_acc00 += re_Apm0*re_Apm0 + im_Apm0*im_Apm0;
            re_acc01 += re_Apm0*re_Apm1 + im_Apm0*im_Apm1;
            re_acc11 += re_Apm1*re_Apm1 + im_Apm1*im_Apm1;
            im_acc00 += re_Apm0*im_Apm0 - im_Apm0*re_Apm0;
            im_acc01 += re_Apm0*im_Apm1 - im_Apm0*re_Apm1;
            im_acc11 += re_Apm1*im_Apm1 - im_Apm1*re_Apm1;
        }
        Rf[((m + 0) * 10 + (m + 0)) * 2 + 0] = re_acc00;
        Rf[((m + 0) * 10 + (m + 0)) * 2 + 1] = im_acc00;
        Rf[((m + 0) * 10 + (m + 1)) * 2 + 0] = re_acc01;
        Rf[((m + 0) * 10 + (m + 1)) * 2 + 1] = im_acc01;
        Rf[((m + 1) * 10 + (m + 0)) * 2 + 0] = re_acc01;
        Rf[((m + 1) * 10 + (m + 0)) * 2 + 1] = -im_acc01;
        Rf[((m + 1) * 10 + (m + 1)) * 2 + 0] = re_acc11;
        Rf[((m + 1) * 10 + (m + 1)) * 2 + 1] = im_acc11;

        /* rest */
        for (n = m + 2; n < 10; n += 2)
        {
            re_acc00 = re_acc11 = re_acc01 = re_acc10 = 0.f;
            im_acc00 = im_acc01 = im_acc10 = im_acc11 = 0.f;

            for (p = 0; p < 10; p++)
            {
                re_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 0];
                re_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 0];
                re_Apn0 = Af[(p * 10 + (n + 0)) * 2 + 0];
                re_Apn1 = Af[(p * 10 + (n + 1)) * 2 + 0];
                im_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 1];
                im_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 1];
                im_Apn0 = Af[(p * 10 + (n + 0)) * 2 + 1];
                im_Apn1 = Af[(p * 10 + (n + 1)) * 2 + 1];

                re_acc00 += re_Apm0*re_Apn0 + im_Apm0*im_Apn0;
                re_acc01 += re_Apm0*re_Apn1 + im_Apm0*im_Apn1;
                re_acc10 += re_Apm1*re_Apn0 + im_Apm1*im_Apn0;
                re_acc11 += re_Apm1*re_Apn1 + im_Apm1*im_Apn1;
                im_acc00 += re_Apm0*im_Apn0 - im_Apm0*re_Apn0;
                im_acc01 += re_Apm0*im_Apn1 - im_Apm0*re_Apn1;
                im_acc10 += re_Apm1*im_Apn0 - im_Apm1*re_Apn0;
                im_acc11 += re_Apm1*im_Apn1 - im_Apm1*re_Apn1;
            }
            Rf[((m + 0) * 10 + (n + 0)) * 2 + 0] = re_acc00;
            Rf[((m + 0) * 10 + (n + 0)) * 2 + 1] = im_acc00;
            Rf[((m + 0) * 10 + (n + 1)) * 2 + 0] = re_acc01;
            Rf[((m + 0) * 10 + (n + 1)) * 2 + 1] = im_acc01;
            Rf[((m + 1) * 10 + (n + 0)) * 2 + 0] = re_acc10;
            Rf[((m + 1) * 10 + (n + 0)) * 2 + 1] = im_acc10;
            Rf[((m + 1) * 10 + (n + 1)) * 2 + 0] = re_acc11;
            Rf[((m + 1) * 10 + (n + 1)) * 2 + 1] = im_acc11;

            Rf[((m + 0) + (n + 0) * 10) * 2 + 0] = re_acc00;
            Rf[((m + 0) + (n + 0) * 10) * 2 + 1] = -im_acc00;
            Rf[((m + 0) + (n + 1) * 10) * 2 + 0] = re_acc01;
            Rf[((m + 0) + (n + 1) * 10) * 2 + 1] = -im_acc01;
            Rf[((m + 1) + (n + 0) * 10) * 2 + 0] = re_acc10;
            Rf[((m + 1) + (n + 0) * 10) * 2 + 1] = -im_acc10;
            Rf[((m + 1) + (n + 1) * 10) * 2 + 0] = re_acc11;
            Rf[((m + 1) + (n + 1) * 10) * 2 + 1] = -im_acc11;

        }
    }
#else
    const float32_t* restrict Af = (const float32_t*)A;
    float32_t* restrict Rf = (float32_t*)R;
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    xtfloat * restrict pRv;
    xtfloat * restrict pRh;
    const xtfloat * restrict pApm;
    const xtfloat * restrict pApn;
    xtfloat R_re, R_im;
    xtfloat Apm_re, Apm_im;
    xtfloat Apn_re, Apn_im;

    for (m = 0; m < 10; m++)
    {
        R_re = sigma2;
        pApm = (xtfloat*)(Af + (0 * 10 + m) * 2);
        pRv = (xtfloat*)(Rf + (m * 10 + m) * 2);
        pRh = (xtfloat*)(Rf + (m + (m + 1) * 10) * 2);

        for (p = 0; p < 10; p++)
        {
            XT_LSIP(Apm_re, pApm, SZ_F32);
            XT_LSIP(Apm_im, pApm, (10 * 2 - 1) * SZ_F32);
            XT_MADD_S(R_re, Apm_re, Apm_re);
            XT_MADD_S(R_re, Apm_im, Apm_im);
        }
        XT_SSIP(R_re, pRv, SZ_F32);
        XT_SSIP(XT_CONST_S(0), pRv, SZ_F32);

        for (n = m + 1; n < 10; n++)
        {
            pApm = (xtfloat*)(Af + (0 * 10 + m) * 2);
            pApn = (xtfloat*)(Af + (0 * 10 + n) * 2);

            R_re = R_im = XT_CONST_S(0);
            for (p = 0; p < 10; p++)
            {
                XT_LSIP(Apm_re, pApm, SZ_F32);
                XT_LSIP(Apm_im, pApm, (10 * 2 - 1) * SZ_F32);
                XT_LSIP(Apn_re, pApn, SZ_F32);
                XT_LSIP(Apn_im, pApn, (10 * 2 - 1) * SZ_F32);
                XT_MADD_S(R_re, Apm_re, Apn_re);
                XT_MADD_S(R_re, Apm_im, Apn_im);
                XT_MADD_S(R_im, Apm_re, Apn_im);
                XT_MSUB_S(R_im, Apm_im, Apn_re);
            }
            XT_SSIP(R_re, pRv, SZ_F32);
            XT_SSIP(R_im, pRv, SZ_F32);
            XT_SSIP(R_re, pRh, SZ_F32);
            XT_SSIP(XT_NEG_S(R_im), pRh, (10 * 2 - 1) * SZ_F32);
        }
    }

#endif
}

size_t  cmatcholpreprocess10x10f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, cmatcholpreprocess10x10f, (void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t  sigma2))

	size_t  cmatcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#endif
