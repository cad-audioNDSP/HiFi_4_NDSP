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

#if (HAVE_VFPU || HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
#define CF32_WIDTH (HIFI_SIMD_WIDTH/SZ_CF32)
/*-------------------------------------------------------------------------
Cholesky MMSE Solver
Compute the MMSE solution for a system of linear equations A*x=B, where A is
an MxN real (complex) matrix with M>=N and rank(A)==N, x is an Nx1 vector of
unknowns, and B is an MxP right hand side vector. This task is accomplished
in 3 steps:
-   Cholesky decomposition is applied to the matrix of normal equations
system, which results in an upper triangular matrix R[NxN] with real and
positive numbers on the main diagonal
-   Forward substitution step: solve R'*y=A'*B for Nx1 vector y.
-   Backward substitution step: solve the system R*x=y for the Nx1 vector of
unknowns x.
For a single MxN matrix A, these 3 steps may be done simultaneously for P
variants of Mx1 right hand side column vectors b gathered into an MxP input
matrix B. MMSE solution is computed independently for each of P columns,
with resulting column vectors forming the solution matrix x of size NxP.
Fixed point API requires explicit setting of fixed point representation of 
input/output matrices as well as for internal temporary matrices such as R 
(Cholesky decomposition) and y (decision of R'y=(A'*B))

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
sigma2      Regularization term. For fixed point, the representation 
            should be Q(2*qA-30)
qRA         qR-qA; difference between fixed point representations of R
            and A matrices (for the fixed point API only). Should be 
            equal or less than 0 (typically -2).
qYBRA       qY-qB+qR-qA, combination of fixed point representations of 
            matrices y,B,R and A (for the fixed point API only)
qXYR        combination of fixed point representation (matrices R, x and y) 
            qXYR=qX-qY+qR
A[M*N]      matrix A. . For fixed point, the representation should be Q(qA)
B[M*P]      Original right-side matrix B. For fixed point, the representation 
            should be Q(qB)
Output:
x[N*P]      Decision matrix x. For fixed point, the representation 
            is Q(qX)
Temporary:
pScr        Scratch data

N = M = 4, 6, 8, 10
P = 1

Restrictions:
All matrices should not overlap and be aligned on 8-bytes
boundary
---------------------------------------------------------------------------*/
void  cmatcholmmsesolver4x4f(void * pScr,
    complex_float * x,
    const complex_float * A,
    const complex_float * B,
    const float32_t       sigma2)
{
    complex_float * restrict D; int SD;
    complex_float * restrict R; int SR;
    complex_float * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(B);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(B, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 4;
    SR = (((4 + 1)*4) >> 1);
    SY = 4*1;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = SD * SZ_CF32;
    SR = SR * SZ_CF32;
    SY = SY * SZ_CF32;
    R = (complex_float *)pScr;
    D = (complex_float *)(((uintptr_t)R) + SR);
    y = (complex_float *)(((uintptr_t)D) + SD);
    pScr = (complex_float *)(((uintptr_t)y) + SY);
    cmatcholdecomp4x4f(pScr, R, D, A, sigma2);
    cmatcholfwdsubst4x4f(pScr, y, R, D, A, B);
    cmatcholbkwsubst4x4f(pScr, x, R, D, y);
}

size_t  cmatcholmmsesolver4x4f_getScratchSize()
{
    size_t s_dc, s_fw, s_bk, s_max;
    int SD, SR, SY;

    SD = 4;
    SR = (((4 + 1)*4) >> 1);
    SY = 4*1;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = SD * SZ_CF32;
    SR = SR * SZ_CF32;
    SY = SY * SZ_CF32;

    s_dc = cmatcholdecomp4x4f_getScratchSize();
    s_fw = cmatcholfwdsubst4x4f_getScratchSize();
    s_bk = cmatcholbkwsubst4x4f_getScratchSize();
    s_max = s_dc > s_fw ? s_dc : s_fw;
    s_max = s_max > s_bk ? s_max : s_bk;
    return SD + SR + SY + s_max;
}
#else
DISCARD_FUN(void, cmatcholmmsesolver4x4f, (void * pScr,
	complex_float * x,
	const complex_float * A,
	const complex_float * B,
	const float32_t sigma2))

size_t  cmatcholmmsesolver4x4f_getScratchSize()
{
	return 0;
}
#endif
