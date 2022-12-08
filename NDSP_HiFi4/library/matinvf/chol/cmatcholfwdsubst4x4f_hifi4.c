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
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*
compute matrix product Z[Nx1]=A[MxN]'*B[Mx1]
Input:
A[SA]    complex matrix MxN
B[SB]    complex matrix Mx1
Output:
Z[N*1]   complex matrix Nx1
*/
void cplxcomputeAB4f(complex_float * Z,
	const complex_float * A,
	const complex_float * B)
{
	xtfloatx2 * restrict pZ = (xtfloatx2 *)Z;
	const xtfloatx2 * restrict pA = (const xtfloatx2 *)A;
	const xtfloatx2 * restrict pB = (const xtfloatx2 *)B;

	xtfloatx2 A0, A1, A2, A3;
	xtfloatx2 B0;
	xtfloatx2 Z0, Z1, Z2, Z3;
	xtfloatx2 Z0t, Z1t, Z2t, Z3t;
	
	Z0 = Z1 = Z2 = Z3 = XT_CONST_S(0);
	Z0t = Z1t = Z2t = Z3t = XT_CONST_S(0);
	// unrolled
	{
		XT_LSX2IP(A0, pA, SZ_CF32);
		XT_LSX2IP(A1, pA, SZ_CF32);
		XT_LSX2IP(A2, pA, SZ_CF32);
		XT_LSX2IP(A3, pA, SZ_CF32);
		XT_LSX2IP(B0, pB, SZ_CF32);
		XT_MADDMUX_S(Z0, A0, B0, 0);
		XT_MADDMUX_S(Z0t, A0, B0, 3);
		XT_MADDMUX_S(Z1, A1, B0, 0);
		XT_MADDMUX_S(Z1t, A1, B0, 3);
		XT_MADDMUX_S(Z2, A2, B0, 0);
		XT_MADDMUX_S(Z2t, A2, B0, 3);
		XT_MADDMUX_S(Z3, A3, B0, 0);
		XT_MADDMUX_S(Z3t, A3, B0, 3);

		XT_LSX2IP(A0, pA, SZ_CF32);
		XT_LSX2IP(A1, pA, SZ_CF32);
		XT_LSX2IP(A2, pA, SZ_CF32);
		XT_LSX2IP(A3, pA, SZ_CF32);
		XT_LSX2IP(B0, pB, SZ_CF32);
		XT_MADDMUX_S(Z0, A0, B0, 0);
		XT_MADDMUX_S(Z0t, A0, B0, 3);
		XT_MADDMUX_S(Z1, A1, B0, 0);
		XT_MADDMUX_S(Z1t, A1, B0, 3);
		XT_MADDMUX_S(Z2, A2, B0, 0);
		XT_MADDMUX_S(Z2t, A2, B0, 3);
		XT_MADDMUX_S(Z3, A3, B0, 0);
		XT_MADDMUX_S(Z3t, A3, B0, 3);

		XT_LSX2IP(A0, pA, SZ_CF32);
		XT_LSX2IP(A1, pA, SZ_CF32);
		XT_LSX2IP(A2, pA, SZ_CF32);
		XT_LSX2IP(A3, pA, SZ_CF32);
		XT_LSX2IP(B0, pB, SZ_CF32);
		XT_MADDMUX_S(Z0, A0, B0, 0);
		XT_MADDMUX_S(Z0t, A0, B0, 3);
		XT_MADDMUX_S(Z1, A1, B0, 0);
		XT_MADDMUX_S(Z1t, A1, B0, 3);
		XT_MADDMUX_S(Z2, A2, B0, 0);
		XT_MADDMUX_S(Z2t, A2, B0, 3);
		XT_MADDMUX_S(Z3, A3, B0, 0);
		XT_MADDMUX_S(Z3t, A3, B0, 3);

		XT_LSX2IP(A0, pA, SZ_CF32);
		XT_LSX2IP(A1, pA, SZ_CF32);
		XT_LSX2IP(A2, pA, SZ_CF32);
		XT_LSX2IP(A3, pA, SZ_CF32);
		XT_LSX2IP(B0, pB, SZ_CF32);
		XT_MADDMUX_S(Z0, A0, B0, 0);
		XT_MADDMUX_S(Z0t, A0, B0, 3);
		XT_MADDMUX_S(Z1, A1, B0, 0);
		XT_MADDMUX_S(Z1t, A1, B0, 3);
		XT_MADDMUX_S(Z2, A2, B0, 0);
		XT_MADDMUX_S(Z2t, A2, B0, 3);
		XT_MADDMUX_S(Z3, A3, B0, 0);
		XT_MADDMUX_S(Z3t, A3, B0, 3);
	}
	Z0 = XT_ADD_SX2(Z0, Z0t);
	Z1 = XT_ADD_SX2(Z1, Z1t);
	Z2 = XT_ADD_SX2(Z2, Z2t);
	Z3 = XT_ADD_SX2(Z3, Z3t);

	XT_SSX2IP(Z0, pZ, SZ_CF32);
	XT_SSX2IP(Z1, pZ, SZ_CF32);
	XT_SSX2IP(Z2, pZ, SZ_CF32);
	XT_SSX2IP(Z3, pZ, SZ_CF32);
	
}
/*-------------------------------------------------------------------------
Cholesky Forward Substitution for Pseudo-inversion
These functions make forward recursion stage of pseudo-inversion. They use
Cholesky decomposition R[NxN] of original matrices A[MxN]

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
            upper triangular matrix R. For fixed point, representation is Q(qR)
A[M*N]      matrix A. For fixed point, representation is Q(qA)
B[M*P]      original right-side matrix B. For fixed point, representation is Q(qB)
D[N]        reciprocals of main diagonal. NOTE: for the fixed point API,
            these data are stored internally in special format with separate
            mantissa and exponent for better accuracy and dynamic range 
            control. So, even for the real data, they stored as pairs of 2
            integers and packed to the complex_fract32 format
qYBRA       qY-qB+qR-qA, combination of fixed point representations of 
            matrices y,B,R and A (for the fixed point API only)
Output:
y[N*P]      Decision matrix y. For fixed point, representation is Q(qY)
Temporary:
pScr        Scratch memory

N = M = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholfwdsubst4x4f(void * pScr,
	complex_float * y,
	const complex_float * R,
	const complex_float * D,
	const complex_float * A,
	const complex_float * B)
{
	NASSERT(y);
	NASSERT(R);
	NASSERT(D);
	NASSERT(A);
	NASSERT(B);
	NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(B, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

	complex_float * _Z = (complex_float *)pScr;

	// compute A'*B
	cplxcomputeAB4f(_Z, A, B);

	// forward recursion
	cplxcholFwdrec4f(y, R, D, _Z, 1);
}

size_t  cmatcholfwdsubst4x4f_getScratchSize()
{
	return 2 * 4 * sizeof(float32_t);
}
#elif (HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*
compute matrix product Z[Nx1]=A[MxN]'*B[Mx1]
Input:
A[SA]    complex matrix MxN
B[SB]    complex matrix Mx1
Output:
Z[N*1]   complex matrix Nx1
*/
void cplxcomputeAB4f(complex_float * Z,
    const complex_float * A,
    const complex_float * B)
{
    float32_t* restrict _Z = (float32_t*)Z;
    const float32_t* restrict _B = (const float32_t*)B;
    const float32_t* restrict _A = (const float32_t*)A;
    float32_t B_re, B_im;
    int n, m;

    for (n = 0; n<4; n++)
    {
        B_re = B_im = 0;
        for (m = 0; m<4; m++)
        {
            float32_t a_re, a_im, b_re, b_im;
            a_re = _A[2 * n + m * 4 * 2 + 0]; a_im = _A[2 * n + m * 4 * 2 + 1];
            b_re = _B[m * 1 * 2 + 0]; b_im = _B[m * 1 * 2 + 1];
            B_re += (a_re*b_re) + (a_im*b_im);
            B_im += (a_re*b_im) - (a_im*b_re);
        }
        _Z[2 * n + 0] = B_re;
        _Z[2 * n + 1] = B_im;
    }
}
/*-------------------------------------------------------------------------
Cholesky Forward Substitution for Pseudo-inversion
These functions make forward recursion stage of pseudo-inversion. They use
Cholesky decomposition R[NxN] of original matrices A[MxN]

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
upper triangular matrix R. For fixed point, representation is Q(qR)
A[M*N]      matrix A. For fixed point, representation is Q(qA)
B[M*P]      original right-side matrix B. For fixed point, representation is Q(qB)
D[N]        reciprocals of main diagonal. NOTE: for the fixed point API,
these data are stored internally in special format with separate
mantissa and exponent for better accuracy and dynamic range
control. So, even for the real data, they stored as pairs of 2
integers and packed to the complex_fract32 format
qYBRA       qY-qB+qR-qA, combination of fixed point representations of
matrices y,B,R and A (for the fixed point API only)
Output:
y[N*P]      Decision matrix y. For fixed point, representation is Q(qY)
Temporary:
pScr        Scratch memory

N = M = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholfwdsubst4x4f(void * pScr,
    complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * A,
    const complex_float * B)
{
    NASSERT(y);
    NASSERT(R);
    NASSERT(D);
    NASSERT(A);
    NASSERT(B);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(B, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    complex_float * _Z = (complex_float *)pScr;

    // compute A'*B
    cplxcomputeAB4f(_Z, A, B);

    // forward recursion
    cplxcholFwdrec4f(y, R, D, _Z, 1);
}

size_t  cmatcholfwdsubst4x4f_getScratchSize()
{
    return 2 * 4 * sizeof(float32_t);
}

#else
DISCARD_FUN(void, cmatcholfwdsubst4x4f, (void * pScr, 
	complex_float * y,
	const complex_float * R, 
	const complex_float * D,
	const complex_float * A, 
	const complex_float * B))

size_t  cmatcholfwdsubst4x4f_getScratchSize()
{
	return 0;
}
#endif
