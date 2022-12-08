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

/* Load 32bit and replicate to xtfloatx2*/
#define _L32_SX2_XP(a,b,c) \
{ \
ae_int32x2 tmp; \
AE_L32_XP(tmp, castxcc(ae_int32, b), c); \
a = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp); }
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
void  matcholpreprocess10x10f(void* pScr,
    float32_t *R,
    const float32_t * A,
    const float32_t  sigma2)
{
    int m, n;
    xtfloatx2 * restrict pRv0;
    //xtfloatx2 * restrict pRv1;
    xtfloatx2 * restrict pRh0;
    //xtfloatx2 * restrict pRh1;
    const xtfloatx2 * restrict pApm0;
    const xtfloatx2 * restrict pApm1;
    const xtfloatx2 * restrict pApn;
    xtfloatx2 Apm0, Apm1;
    xtfloatx2 Apn01;
    xtfloatx2 Acc0, Acc1, Acc0t, Acc1t;
    xtfloatx2 Acc0tt, Acc1tt;
    xtfloatx2 sigma = sigma2;
    xtfloatx2 sigma_w_zeros, zeros_w_sigma;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    sigma_w_zeros = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(AE_SLAI64(AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(sigma)), 32)));
    zeros_w_sigma = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(AE_SRAI64(AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(sigma)), 32)));

    for (m = 0; m < 10; m += 2)
    {
        Acc0 = sigma_w_zeros;
        Acc1 = zeros_w_sigma;
        pRv0 = pRh0 = (xtfloatx2*)(R + m * (10 + 1));
        pApm0 = pApn = (xtfloatx2*)(A + m);
        pApm1 = (xtfloatx2*)(A + m + 1);
        for (n = m; n < 10; n += 2)
        {
            //unrolled 
            //for (p = 0; p < 5; p++)
            {
                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0, Apm0, Apn01);
                XT_MADD_SX2(Acc1, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                Acc0t = XT_MUL_SX2(Apm0, Apn01);
                Acc1t = XT_MUL_SX2(Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                Acc0tt = XT_MUL_SX2(Apm0, Apn01);
                Acc1tt = XT_MUL_SX2(Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0, Apm0, Apn01);
                XT_MADD_SX2(Acc1, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0t, Apm0, Apn01);
                XT_MADD_SX2(Acc1t, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0tt, Apm0, Apn01);
                XT_MADD_SX2(Acc1tt, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0, Apm0, Apn01);
                XT_MADD_SX2(Acc1, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0t, Apm0, Apn01);
                XT_MADD_SX2(Acc1t, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, 10 * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, 10 * SZ_F32);
                XT_LSX2IP(Apn01, pApn, 10 * SZ_F32);

                XT_MADD_SX2(Acc0tt, Apm0, Apn01);
                XT_MADD_SX2(Acc1tt, Apm1, Apn01);

                _L32_SX2_XP(Apm0, pApm0, -(9 * 10) * SZ_F32);
                _L32_SX2_XP(Apm1, pApm1, -(9 * 10) * SZ_F32);
                XT_LSX2XP(Apn01, pApn, -(9 * 10 - 2) * SZ_F32);

                XT_MADD_SX2(Acc0, Apm0, Apn01);
                XT_MADD_SX2(Acc1, Apm1, Apn01);
            }
            Acc0t = XT_ADD_SX2(Acc0t, Acc0tt);
            Acc0 = XT_ADD_SX2(Acc0, Acc0t);
            Acc1t = XT_ADD_SX2(Acc1t, Acc1tt);
            Acc1 = XT_ADD_SX2(Acc1, Acc1t);

            XT_SSX2IP(Acc0, pRv0, 10 * SZ_F32);
            XT_SSX2XP(Acc1, pRv0, -(10 - 2) * SZ_F32);
            XT_SSX2IP(XT_SEL32_HH_SX2(Acc0, Acc1), pRh0, 10 * SZ_F32);
            XT_SSX2IP(XT_SEL32_LL_SX2(Acc0, Acc1), pRh0, 10 * SZ_F32);

            Acc0 = Acc1 = XT_CONST_S(0);
        }
    }
}
/* scratch allocation functions */
size_t   matcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#elif(HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
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
/* single-matrix API */
void  matcholpreprocess10x10f(void* restrict pScr,
    float32_t * restrict R,
    const float32_t * restrict A,
    const float32_t  sigma2)
{
#if 0
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    float32_t acc;
    float32_t Apm;
    float32_t Apn;
    for (m = 0; m < 10; m++)
    {
        acc = sigma2;
        for (p = 0; p < 10; p++)
        {
            Apm = A[p * 10 + m];

            acc += Apm*Apm;
        }
        R[m * 10 + m] = acc;

        for (n = m + 1; n < 10; n++)
        {
            acc = 0.f;
            for (p = 0; p < 10; p++)
            {
                Apm = A[p * 10 + m];
                Apn = A[p * 10 + n];

                acc += Apm*Apn;
            }
            R[m * 10 + n] = acc;
            R[m + n * 10] = acc;
        }
    }
#else
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    float32_t acc00, acc01, acc10, acc11;
    float32_t Apm0, Apm1;
    float32_t Apn0, Apn1;
    /* processing by squares 2x2 */
    for (m = 0; m < 10; m+=2)
    {
        /* diagonal squares */
        acc00 = acc11 = sigma2;
        acc01 = 0.f;
        for (p = 0; p < 10; p++)
        {
            Apm0 = A[p * 10 + m + 0];
            Apm1 = A[p * 10 + m + 1];
            acc00 += Apm0*Apm0;
            acc01 += Apm0*Apm1;
            acc11 += Apm1*Apm1;
        }
        R[(m + 0) * 10 + (m + 0)] = acc00;
        R[(m + 0) * 10 + (m + 1)] = acc01;
        R[(m + 1) * 10 + (m + 0)] = acc01;
        R[(m + 1) * 10 + (m + 1)] = acc11;

        /* rest */
        for (n = m+2; n < 10; n+=2)
        {
            acc00 = acc11 = 0.f;
            acc01 = acc10 = 0.f;
            for (p = 0; p < 10; p++)
            {
                Apm0 = A[p * 10 + m + 0];
                Apm1 = A[p * 10 + m + 1];
                Apn0 = A[p * 10 + n + 0];
                Apn1 = A[p * 10 + n + 1];
                acc00 += Apm0*Apn0;
                acc01 += Apm0*Apn1;
                acc10 += Apm1*Apn0;
                acc11 += Apm1*Apn1;
            }
            R[(m + 0) * 10 + (n + 0)] = acc00;
            R[(m + 0) * 10 + (n + 1)] = acc01;
            R[(m + 1) * 10 + (n + 0)] = acc10;
            R[(m + 1) * 10 + (n + 1)] = acc11;

            R[(m + 0) + (n + 0) * 10] = acc00;
            R[(m + 0) + (n + 1) * 10] = acc01;
            R[(m + 1) + (n + 0) * 10] = acc10;
            R[(m + 1) + (n + 1) * 10] = acc11;
        }
    }
#endif
}

size_t  matcholpreprocess10x10f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, matcholpreprocess10x10f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

	size_t  matcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#endif
