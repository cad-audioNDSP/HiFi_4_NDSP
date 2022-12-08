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
#define SZ_F32 (sizeof(float32_t))

/* MOVT64 for xtfloatx2 */
#define _MOVT64_SX2(a,b,c)  \
{ \
ae_int64 tmp = AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(a)); \
AE_MOVT64(tmp, AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(b)), c); \
a = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(tmp)); }

/* Load 32bit and replicate to xtfloatx2*/
#define _L32_SX2_IP(a,b,c) \
{ \
ae_int32x2 tmp; \
AE_L32_IP(tmp, castxcc(ae_int32, b), c); \
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
void  matcholpreprocess6x6f(void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2)
{
	const xtfloatx2 * restrict pApn;
	const xtfloatx2 * restrict pApm;
	xtfloatx2 * restrict pR = (xtfloatx2 *)R;
	xtfloatx2 Apm;
	xtfloatx2 Apn0, Apn1, Apn2, Apn3, Apn4, Apn5;
	xtfloatx2 Acc0, Acc1, Acc2, Acc3, Acc4, Acc5;
	xtfloatx2 sigma = sigma2;
	xtfloatx2 sigma_w_zeros;
	xtfloatx2 zeros_w_sigma;
	int p;

	NASSERT(R);
	NASSERT(A);
	NASSERT(sigma2);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

	sigma_w_zeros = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(AE_SLAI64(AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(sigma)), 32)));
	zeros_w_sigma = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(AE_SRAI64(AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(sigma)), 32)));

	/* fully unrolled */
	{
		pApm = (xtfloatx2 *)(A + 0);
		pApn = (xtfloatx2 *)(A);

		Acc0 = sigma_w_zeros;
		Acc1 = Acc2 = Acc3 = Acc4 = Acc5 = XT_CONST_S(0);

		/* unrolled 1 time */
		for (p = 0; p < 3; p++)
		{
			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn0, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn1, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn2, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc0, Apm, Apn0);
			XT_MADD_SX2(Acc1, Apm, Apn1);
			XT_MADD_SX2(Acc2, Apm, Apn2);

			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn3, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn4, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn5, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc3, Apm, Apn3);
			XT_MADD_SX2(Acc4, Apm, Apn4);
			XT_MADD_SX2(Acc5, Apm, Apn5);
		}
		Acc0 = XT_ADD_SX2(Acc0, Acc3);
		Acc1 = XT_ADD_SX2(Acc1, Acc4);
		Acc2 = XT_ADD_SX2(Acc2, Acc5);

		XT_SSX2IP(Acc0, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc1, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc2, pR, 2 * SZ_F32);

		pApm = (xtfloatx2 *)(A + 1);
		pApn = (xtfloatx2 *)(A);

		Acc0 = zeros_w_sigma;
		Acc1 = Acc2 = Acc3 = Acc4 = Acc5 = XT_CONST_S(0);

		/* unrolled 1 time */
		for (p = 0; p < 3; p++)
		{
			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn0, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn1, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn2, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc0, Apm, Apn0);
			XT_MADD_SX2(Acc1, Apm, Apn1);
			XT_MADD_SX2(Acc2, Apm, Apn2);

			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn3, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn4, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn5, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc3, Apm, Apn3);
			XT_MADD_SX2(Acc4, Apm, Apn4);
			XT_MADD_SX2(Acc5, Apm, Apn5);
		}
		Acc0 = XT_ADD_SX2(Acc0, Acc3);
		Acc1 = XT_ADD_SX2(Acc1, Acc4);
		Acc2 = XT_ADD_SX2(Acc2, Acc5);

		XT_SSX2IP(Acc0, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc1, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc2, pR, 2 * SZ_F32);

		pApm = (xtfloatx2 *)(A + 2);
		pApn = (xtfloatx2 *)(A);

		Acc1 = sigma_w_zeros;
		Acc0 = Acc2 = Acc3 = Acc4 = Acc5 = XT_CONST_S(0);

		/* unrolled 1 time */
		for (p = 0; p < 3; p++)
		{
			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn0, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn1, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn2, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc0, Apm, Apn0);
			XT_MADD_SX2(Acc1, Apm, Apn1);
			XT_MADD_SX2(Acc2, Apm, Apn2);

			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn3, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn4, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn5, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc3, Apm, Apn3);
			XT_MADD_SX2(Acc4, Apm, Apn4);
			XT_MADD_SX2(Acc5, Apm, Apn5);
		}
		Acc0 = XT_ADD_SX2(Acc0, Acc3);
		Acc1 = XT_ADD_SX2(Acc1, Acc4);
		Acc2 = XT_ADD_SX2(Acc2, Acc5);

		XT_SSX2IP(Acc0, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc1, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc2, pR, 2 * SZ_F32);

		pApm = (xtfloatx2 *)(A + 3);
		pApn = (xtfloatx2 *)(A);

		Acc1 = zeros_w_sigma;
		Acc0 = Acc2 = Acc3 = Acc4 = Acc5 = XT_CONST_S(0);

		/* unrolled 1 time */
		for (p = 0; p < 3; p++)
		{
			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn0, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn1, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn2, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc0, Apm, Apn0);
			XT_MADD_SX2(Acc1, Apm, Apn1);
			XT_MADD_SX2(Acc2, Apm, Apn2);

			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn3, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn4, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn5, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc3, Apm, Apn3);
			XT_MADD_SX2(Acc4, Apm, Apn4);
			XT_MADD_SX2(Acc5, Apm, Apn5);
		}
		Acc0 = XT_ADD_SX2(Acc0, Acc3);
		Acc1 = XT_ADD_SX2(Acc1, Acc4);
		Acc2 = XT_ADD_SX2(Acc2, Acc5);

		XT_SSX2IP(Acc0, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc1, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc2, pR, 2 * SZ_F32);

		pApm = (xtfloatx2 *)(A + 4);
		pApn = (xtfloatx2 *)(A);

		Acc2 = sigma_w_zeros;
		Acc0 = Acc1 = Acc3 = Acc4 = Acc5 = XT_CONST_S(0);

		/* unrolled 1 time */
		for (p = 0; p < 3; p++)
		{
			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn0, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn1, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn2, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc0, Apm, Apn0);
			XT_MADD_SX2(Acc1, Apm, Apn1);
			XT_MADD_SX2(Acc2, Apm, Apn2);

			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn3, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn4, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn5, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc3, Apm, Apn3);
			XT_MADD_SX2(Acc4, Apm, Apn4);
			XT_MADD_SX2(Acc5, Apm, Apn5);
		}
		Acc0 = XT_ADD_SX2(Acc0, Acc3);
		Acc1 = XT_ADD_SX2(Acc1, Acc4);
		Acc2 = XT_ADD_SX2(Acc2, Acc5);

		XT_SSX2IP(Acc0, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc1, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc2, pR, 2 * SZ_F32);

		pApm = (xtfloatx2 *)(A + 5);
		pApn = (xtfloatx2 *)(A);

		Acc2 = zeros_w_sigma;
		Acc0 = Acc1 = Acc3 = Acc4 = Acc5 = XT_CONST_S(0);

		/* unrolled 1 time */
		for (p = 0; p < 3; p++)
		{
			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn0, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn1, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn2, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc0, Apm, Apn0);
			XT_MADD_SX2(Acc1, Apm, Apn1);
			XT_MADD_SX2(Acc2, Apm, Apn2);

			_L32_SX2_IP(Apm, pApm, 6 * SZ_F32);
			XT_LSX2IP(Apn3, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn4, pApn, 2 * SZ_F32);
			XT_LSX2IP(Apn5, pApn, 2 * SZ_F32);
			XT_MADD_SX2(Acc3, Apm, Apn3);
			XT_MADD_SX2(Acc4, Apm, Apn4);
			XT_MADD_SX2(Acc5, Apm, Apn5);
		}
		Acc0 = XT_ADD_SX2(Acc0, Acc3);
		Acc1 = XT_ADD_SX2(Acc1, Acc4);
		Acc2 = XT_ADD_SX2(Acc2, Acc5);

		XT_SSX2IP(Acc0, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc1, pR, 2 * SZ_F32);
		XT_SSX2IP(Acc2, pR, 2 * SZ_F32);
	}
}

/* scratch allocation functions */
size_t   matcholpreprocess6x6f_getScratchSize()
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
void  matcholpreprocess6x6f(void* pScr,
    float32_t * R,
    const float32_t * A,
    const float32_t  sigma2)
{
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
    for (m = 0; m < 6; m += 2)
    {
        /* diagonal squares */
        acc00 = acc11 = sigma2;
        acc01 = 0.f;
        for (p = 0; p < 6; p++)
        {
            Apm0 = A[p * 6 + m + 0];
            Apm1 = A[p * 6 + m + 1];
            acc00 += Apm0*Apm0;
            acc01 += Apm0*Apm1;
            acc11 += Apm1*Apm1;
        }
        R[(m + 0) * 6 + (m + 0)] = acc00;
        R[(m + 0) * 6 + (m + 1)] = acc01;
        R[(m + 1) * 6 + (m + 0)] = acc01;
        R[(m + 1) * 6 + (m + 1)] = acc11;

        /* rest */
        for (n = m + 2; n < 6; n += 2)
        {
            acc00 = acc11 = 0.f;
            acc01 = acc10 = 0.f;
            for (p = 0; p < 6; p++)
            {
                Apm0 = A[p * 6 + m + 0];
                Apm1 = A[p * 6 + m + 1];
                Apn0 = A[p * 6 + n + 0];
                Apn1 = A[p * 6 + n + 1];
                acc00 += Apm0*Apn0;
                acc01 += Apm0*Apn1;
                acc10 += Apm1*Apn0;
                acc11 += Apm1*Apn1;
            }
            R[(m + 0) * 6 + (n + 0)] = acc00;
            R[(m + 0) * 6 + (n + 1)] = acc01;
            R[(m + 1) * 6 + (n + 0)] = acc10;
            R[(m + 1) * 6 + (n + 1)] = acc11;

            R[(m + 0) + (n + 0) * 6] = acc00;
            R[(m + 0) + (n + 1) * 6] = acc01;
            R[(m + 1) + (n + 0) * 6] = acc10;
            R[(m + 1) + (n + 1) * 6] = acc11;
        }
    }
}

size_t  matcholpreprocess6x6f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, matcholpreprocess6x6f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

	size_t  matcholpreprocess6x6f_getScratchSize()
{
	return 0;
}
#endif
