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
Cholesky Backward Substitution for Pseudo-inversion
These functions make backward substitution stage of pseudo-inversion. They
use Cholesky decomposition of original matrices and results of forward
substitution.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
            Cholesky upper-triangle matrix R. For fixed point API, the 
            representation is Q(qR)
y[N*P]      Results of forward substitution stage. For fixed point API, 
            the representation is Q(qY)
D[N]        sequence of reciprocals of main diagonal R. NOTE: for the 
            fixed point API, these data are stored internally in special 
            format with separate mantissa and exponent for better accuracy 
            and dynamic range control. So, even for the real data, they 
            stored as pairs of 2 integers and packed to the complex_fract32 
            format
qXYR        combination of fixed point representation (matrices R, x and y) 
            qXYR=qX-qY+qR (for fixed point API only)
Output:
x[N*P]      Decision matrix x. For fixed point API, the representation is Q(qX)

N = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholbkwsubst10x10f(void * pScr,
    complex_float * x,
    const complex_float * R,
    const complex_float * D,
    const complex_float * y)
{
    NASSERT(x);
    NASSERT(R);
    NASSERT(D);
    NASSERT(y);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    //NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    const xtfloatx2 * restrict pR9 = (const xtfloatx2*)(R + (10 * 9) / 2 + 0); //N*(N-1)/2
    const xtfloatx2 * restrict pR8 = (const xtfloatx2*)(R + (10 * 9) / 2 + 1);
    const xtfloatx2 * restrict pR7 = (const xtfloatx2*)(R + (10 * 9) / 2 + 2);
    const xtfloatx2 * restrict pR6 = (const xtfloatx2*)(R + (10 * 9) / 2 + 3);
    const xtfloatx2 * restrict pR5 = (const xtfloatx2*)(R + (10 * 9) / 2 + 4);
    const xtfloatx2 * restrict pR4 = (const xtfloatx2*)(R + (10 * 9) / 2 + 5);
    const xtfloatx2 * restrict pR3 = (const xtfloatx2*)(R + (10 * 9) / 2 + 6);
    const xtfloatx2 * restrict pR2 = (const xtfloatx2*)(R + (10 * 9) / 2 + 7);
    const xtfloatx2 * restrict pR1 = (const xtfloatx2*)(R + (10 * 9) / 2 + 8);
    const xtfloatx2 * restrict pY = (const xtfloatx2*)(y + 9);
    const xtfloatx2 * restrict pD = (const xtfloatx2*)(D + 9);
    xtfloatx2 * restrict pX = (xtfloatx2*)(x + 9);
    xtfloatx2 D0;
    xtfloatx2 R1, R2, R3, R4, R5, R6, R7, R8, R9;
    xtfloatx2 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, Xt;

    // X0
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X0, pY, -SZ_CF32);
    X0 = XT_MUL_SX2(X0, D0);
    XT_SSX2XP(X0, pX, -SZ_CF32);

    // X1
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X1, pY, -SZ_CF32);
    XT_LSX2IP(R1, pR1, 0);
    XT_MADDMUX_S(X1, X0, R1, 2);
    XT_MADDMUX_S(X1, X0, R1, 3);
    X1 = XT_MUL_SX2(X1, D0);
    XT_SSX2XP(X1, pX, -SZ_CF32);

    // X2
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X2, pY, -SZ_CF32);
    XT_LSX2XP(R2, pR2, -SZ_CF32 * 9); // N-1
    XT_MADDMUX_S(X2, X0, R2, 2);
    XT_MADDMUX_S(Xt, X0, R2, 3);
    XT_LSX2IP(R2, pR2, 0);
    XT_MADDMUX_S(X2, X1, R2, 2);
    XT_MADDMUX_S(Xt, X1, R2, 3);
    X2 = XT_ADD_SX2(X2, Xt);
    X2 = XT_MUL_SX2(X2, D0);
    XT_SSX2XP(X2, pX, -SZ_CF32);

    // X3
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X3, pY, -SZ_CF32);
    XT_LSX2XP(R3, pR3, -SZ_CF32 * 9);
    XT_MADDMUX_S(X3, X0, R3, 2);
    XT_MADDMUX_S(Xt, X0, R3, 3);
    XT_LSX2XP(R3, pR3, -SZ_CF32 * 8);
    XT_MADDMUX_S(X3, X1, R3, 2);
    XT_MADDMUX_S(Xt, X1, R3, 3);
    XT_LSX2IP(R3, pR3, 0);
    XT_MADDMUX_S(X3, X2, R3, 2);
    XT_MADDMUX_S(Xt, X2, R3, 3);
    X3 = XT_ADD_SX2(X3, Xt);
    X3 = XT_MUL_SX2(X3, D0);
    XT_SSX2XP(X3, pX, -SZ_CF32);

    // X4
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X4, pY, -SZ_CF32);
    XT_LSX2XP(R4, pR4, -SZ_CF32 * 9);
    XT_MADDMUX_S(X4, X0, R4, 2);
    XT_MADDMUX_S(Xt, X0, R4, 3);
    XT_LSX2XP(R4, pR4, -SZ_CF32 * 8);
    XT_MADDMUX_S(X4, X1, R4, 2);
    XT_MADDMUX_S(Xt, X1, R4, 3);
    XT_LSX2XP(R4, pR4, -SZ_CF32 * 7);
    XT_MADDMUX_S(X4, X2, R4, 2);
    XT_MADDMUX_S(Xt, X2, R4, 3);
    XT_LSX2IP(R4, pR4, 0);
    XT_MADDMUX_S(X4, X3, R4, 2);
    XT_MADDMUX_S(Xt, X3, R4, 3);
    X4 = XT_ADD_SX2(X4, Xt);
    X4 = XT_MUL_SX2(X4, D0);
    XT_SSX2XP(X4, pX, -SZ_CF32);

    // X5
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X5, pY, -SZ_CF32);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 9);
    XT_MADDMUX_S(X5, X0, R5, 2);
    XT_MADDMUX_S(Xt, X0, R5, 3);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 8);
    XT_MADDMUX_S(X5, X1, R5, 2);
    XT_MADDMUX_S(Xt, X1, R5, 3);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 7);
    XT_MADDMUX_S(X5, X2, R5, 2);
    XT_MADDMUX_S(Xt, X2, R5, 3);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 6);
    XT_MADDMUX_S(X5, X3, R5, 2);
    XT_MADDMUX_S(Xt, X3, R5, 3);
    XT_LSX2IP(R5, pR5, 0);
    XT_MADDMUX_S(X5, X4, R5, 2);
    XT_MADDMUX_S(Xt, X4, R5, 3);
    X5 = XT_ADD_SX2(X5, Xt);
    X5 = XT_MUL_SX2(X5, D0);
    XT_SSX2XP(X5, pX, -SZ_CF32);

    // X6
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X6, pY, -SZ_CF32);
    XT_LSX2XP(R6, pR6, -SZ_CF32 * 9);
    XT_MADDMUX_S(X6, X0, R6, 2);
    XT_MADDMUX_S(Xt, X0, R6, 3);
    XT_LSX2XP(R6, pR6, -SZ_CF32 * 8);
    XT_MADDMUX_S(X6, X1, R6, 2);
    XT_MADDMUX_S(Xt, X1, R6, 3);
    XT_LSX2XP(R6, pR6, -SZ_CF32 * 7);
    XT_MADDMUX_S(X6, X2, R6, 2);
    XT_MADDMUX_S(Xt, X2, R6, 3);
    XT_LSX2XP(R6, pR6, -SZ_CF32 * 6);
    XT_MADDMUX_S(X6, X3, R6, 2);
    XT_MADDMUX_S(Xt, X3, R6, 3);
    XT_LSX2XP(R6, pR6, -SZ_CF32 * 5);
    XT_MADDMUX_S(X6, X4, R6, 2);
    XT_MADDMUX_S(Xt, X4, R6, 3);
    XT_LSX2IP(R6, pR6, 0);
    XT_MADDMUX_S(X6, X5, R6, 2);
    XT_MADDMUX_S(Xt, X5, R6, 3);
    X6 = XT_ADD_SX2(X6, Xt);
    X6 = XT_MUL_SX2(X6, D0);
    XT_SSX2XP(X6, pX, -SZ_CF32);

    // X7
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X7, pY, -SZ_CF32);
    XT_LSX2XP(R7, pR7, -SZ_CF32 * 9);
    XT_MADDMUX_S(X7, X0, R7, 2);
    XT_MADDMUX_S(Xt, X0, R7, 3);
    XT_LSX2XP(R7, pR7, -SZ_CF32 * 8);
    XT_MADDMUX_S(X7, X1, R7, 2);
    XT_MADDMUX_S(Xt, X1, R7, 3);
    XT_LSX2XP(R7, pR7, -SZ_CF32 * 7);
    XT_MADDMUX_S(X7, X2, R7, 2);
    XT_MADDMUX_S(Xt, X2, R7, 3);
    XT_LSX2XP(R7, pR7, -SZ_CF32 * 6);
    XT_MADDMUX_S(X7, X3, R7, 2);
    XT_MADDMUX_S(Xt, X3, R7, 3);
    XT_LSX2XP(R7, pR7, -SZ_CF32 * 5);
    XT_MADDMUX_S(X7, X4, R7, 2);
    XT_MADDMUX_S(Xt, X4, R7, 3);
    XT_LSX2XP(R7, pR7, -SZ_CF32 * 4);
    XT_MADDMUX_S(X7, X5, R7, 2);
    XT_MADDMUX_S(Xt, X5, R7, 3);
    XT_LSX2IP(R7, pR7, 0);
    XT_MADDMUX_S(X7, X6, R7, 2);
    XT_MADDMUX_S(Xt, X6, R7, 3);
    X7 = XT_ADD_SX2(X7, Xt);
    X7 = XT_MUL_SX2(X7, D0);
    XT_SSX2XP(X7, pX, -SZ_CF32);

    // X8
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X8, pY, -SZ_CF32);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 9);
    XT_MADDMUX_S(X8, X0, R8, 2);
    XT_MADDMUX_S(Xt, X0, R8, 3);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 8);
    XT_MADDMUX_S(X8, X1, R8, 2);
    XT_MADDMUX_S(Xt, X1, R8, 3);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 7);
    XT_MADDMUX_S(X8, X2, R8, 2);
    XT_MADDMUX_S(Xt, X2, R8, 3);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 6);
    XT_MADDMUX_S(X8, X3, R8, 2);
    XT_MADDMUX_S(Xt, X3, R8, 3);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 5);
    XT_MADDMUX_S(X8, X4, R8, 2);
    XT_MADDMUX_S(Xt, X4, R8, 3);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 4);
    XT_MADDMUX_S(X8, X5, R8, 2);
    XT_MADDMUX_S(Xt, X5, R8, 3);
    XT_LSX2XP(R8, pR8, -SZ_CF32 * 3);
    XT_MADDMUX_S(X8, X6, R8, 2);
    XT_MADDMUX_S(Xt, X6, R8, 3);
    XT_LSX2IP(R8, pR8, 0);
    XT_MADDMUX_S(X8, X7, R8, 2);
    XT_MADDMUX_S(Xt, X7, R8, 3);
    X8 = XT_ADD_SX2(X8, Xt);
    X8 = XT_MUL_SX2(X8, D0);
    XT_SSX2XP(X8, pX, -SZ_CF32);

    // X9
    Xt = XT_CONST_S(0);
    XT_LSX2XP(D0, pD, -SZ_CF32);
    XT_LSX2XP(X9, pY, -SZ_CF32);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 9);
    XT_MADDMUX_S(X9, X0, R9, 2);
    XT_MADDMUX_S(Xt, X0, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 8);
    XT_MADDMUX_S(X9, X1, R9, 2);
    XT_MADDMUX_S(Xt, X1, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 7);
    XT_MADDMUX_S(X9, X2, R9, 2);
    XT_MADDMUX_S(Xt, X2, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 6);
    XT_MADDMUX_S(X9, X3, R9, 2);
    XT_MADDMUX_S(Xt, X3, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 5);
    XT_MADDMUX_S(X9, X4, R9, 2);
    XT_MADDMUX_S(Xt, X4, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 4);
    XT_MADDMUX_S(X9, X5, R9, 2);
    XT_MADDMUX_S(Xt, X5, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 3);
    XT_MADDMUX_S(X9, X6, R9, 2);
    XT_MADDMUX_S(Xt, X6, R9, 3);
    XT_LSX2XP(R9, pR9, -SZ_CF32 * 2);
    XT_MADDMUX_S(X9, X7, R9, 2);
    XT_MADDMUX_S(Xt, X7, R9, 3);
    XT_LSX2IP(R9, pR9, 0);
    XT_MADDMUX_S(X9, X8, R9, 2);
    XT_MADDMUX_S(Xt, X8, R9, 3);
    X9 = XT_ADD_SX2(X9, Xt);
    X9 = XT_MUL_SX2(X9, D0);
    XT_SSX2XP(X9, pX, -SZ_CF32);

}

size_t  cmatcholbkwsubst10x10f_getScratchSize()
{
    return 0;
}
#elif (HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*-------------------------------------------------------------------------
Cholesky Backward Substitution for Pseudo-inversion
These functions make backward substitution stage of pseudo-inversion. They
use Cholesky decomposition of original matrices and results of forward
substitution.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
Cholesky upper-triangle matrix R. For fixed point API, the
representation is Q(qR)
y[N*P]      Results of forward substitution stage. For fixed point API,
the representation is Q(qY)
D[N]        sequence of reciprocals of main diagonal R. NOTE: for the
fixed point API, these data are stored internally in special
format with separate mantissa and exponent for better accuracy
and dynamic range control. So, even for the real data, they
stored as pairs of 2 integers and packed to the complex_fract32
format
qXYR        combination of fixed point representation (matrices R, x and y)
qXYR=qX-qY+qR (for fixed point API only)
Output:
x[N*P]      Decision matrix x. For fixed point API, the representation is Q(qX)

N = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholbkwsubst10x10f(void * pScr,
    complex_float * x,
    const complex_float * R,
    const complex_float * D,
    const complex_float * y)
{
    NASSERT(x);
    NASSERT(R);
    NASSERT(D);
    NASSERT(y);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    //NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    const xtfloat * restrict pR9 = (const xtfloat*)(R + (10 * 9) / 2 + 0) + 1; //N*(N-1)/2
    const xtfloat * restrict pR8 = (const xtfloat*)(R + (10 * 9) / 2 + 1) + 1;
    const xtfloat * restrict pR7 = (const xtfloat*)(R + (10 * 9) / 2 + 2) + 1;
    const xtfloat * restrict pR6 = (const xtfloat*)(R + (10 * 9) / 2 + 3) + 1;
    const xtfloat * restrict pR5 = (const xtfloat*)(R + (10 * 9) / 2 + 4) + 1;
    const xtfloat * restrict pR4 = (const xtfloat*)(R + (10 * 9) / 2 + 5) + 1;
    const xtfloat * restrict pR3 = (const xtfloat*)(R + (10 * 9) / 2 + 6) + 1;
    const xtfloat * restrict pR2 = (const xtfloat*)(R + (10 * 9) / 2 + 7) + 1;
    const xtfloat * restrict pR1 = (const xtfloat*)(R + (10 * 9) / 2 + 8) + 1;
    const xtfloat * restrict pY = (const xtfloat*)(y + 9) + 1;
    const xtfloat * restrict pD = (const xtfloat*)(D + 9) + 1;
    xtfloat * restrict pX = (xtfloat*)(x + 9) + 1;
    xtfloat D0_re, D0_im;
    xtfloat R1_re, R2_re, R3_re, R4_re, R5_re, R6_re, R7_re, R8_re, R9_re;
    xtfloat R1_im, R2_im, R3_im, R4_im, R5_im, R6_im, R7_im, R8_im, R9_im;
    xtfloat X0_re, X1_re, X2_re, X3_re, X4_re, X5_re, X6_re, X7_re, X8_re, X9_re;
    xtfloat X0_im, X1_im, X2_im, X3_im, X4_im, X5_im, X6_im, X7_im, X8_im, X9_im;

    // X0
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X0_im, pY, -SZ_F32);
    XT_LSXP(X0_re, pY, -SZ_F32);
    X0_re = XT_MUL_S(X0_re, D0_re);
    X0_im = XT_MUL_S(X0_im, D0_im);
    XT_SSXP(X0_im, pX, -SZ_F32);
    XT_SSXP(X0_re, pX, -SZ_F32);

    // X1
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X1_im, pY, -SZ_F32);
    XT_LSXP(X1_re, pY, -SZ_F32);
    XT_LSXP(R1_im, pR1, -SZ_F32);
    XT_LSIP(R1_re, pR1, 0);
    XT_MSUB_S(X1_re, X0_re, R1_re);
    XT_MADD_S(X1_re, X0_im, R1_im);
    XT_MSUB_S(X1_im, X0_re, R1_im);
    XT_MSUB_S(X1_im, X0_im, R1_re);
    X1_re = XT_MUL_S(X1_re, D0_re);
    X1_im = XT_MUL_S(X1_im, D0_im);
    XT_SSXP(X1_im, pX, -SZ_F32);
    XT_SSXP(X1_re, pX, -SZ_F32);

    // X2
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X2_im, pY, -SZ_F32);
    XT_LSXP(X2_re, pY, -SZ_F32);
    XT_LSXP(R2_im, pR2, -SZ_F32);
    XT_LSXP(R2_re, pR2, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X2_re, X0_re, R2_re);
    XT_MADD_S(X2_re, X0_im, R2_im);
    XT_MSUB_S(X2_im, X0_re, R2_im);
    XT_MSUB_S(X2_im, X0_im, R2_re);
    XT_LSXP(R2_im, pR2, -SZ_F32);
    XT_LSIP(R2_re, pR2, 0);
    XT_MSUB_S(X2_re, X1_re, R2_re);
    XT_MADD_S(X2_re, X1_im, R2_im);
    XT_MSUB_S(X2_im, X1_re, R2_im);
    XT_MSUB_S(X2_im, X1_im, R2_re);
    X2_re = XT_MUL_S(X2_re, D0_re);
    X2_im = XT_MUL_S(X2_im, D0_im);
    XT_SSXP(X2_im, pX, -SZ_F32);
    XT_SSXP(X2_re, pX, -SZ_F32);

    // X3
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X3_im, pY, -SZ_F32);
    XT_LSXP(X3_re, pY, -SZ_F32);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSXP(R3_re, pR3, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X3_re, X0_re, R3_re);
    XT_MADD_S(X3_re, X0_im, R3_im);
    XT_MSUB_S(X3_im, X0_re, R3_im);
    XT_MSUB_S(X3_im, X0_im, R3_re);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSXP(R3_re, pR3, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X3_re, X1_re, R3_re);
    XT_MADD_S(X3_re, X1_im, R3_im);
    XT_MSUB_S(X3_im, X1_re, R3_im);
    XT_MSUB_S(X3_im, X1_im, R3_re);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSIP(R3_re, pR3, 0);
    XT_MSUB_S(X3_re, X2_re, R3_re);
    XT_MADD_S(X3_re, X2_im, R3_im);
    XT_MSUB_S(X3_im, X2_re, R3_im);
    XT_MSUB_S(X3_im, X2_im, R3_re);
    X3_re = XT_MUL_S(X3_re, D0_re);
    X3_im = XT_MUL_S(X3_im, D0_im);
    XT_SSXP(X3_im, pX, -SZ_F32);
    XT_SSXP(X3_re, pX, -SZ_F32);

    // X4
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X4_im, pY, -SZ_F32);
    XT_LSXP(X4_re, pY, -SZ_F32);
    XT_LSXP(R4_im, pR4, -SZ_F32);
    XT_LSXP(R4_re, pR4, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X4_re, X0_re, R4_re);
    XT_MADD_S(X4_re, X0_im, R4_im);
    XT_MSUB_S(X4_im, X0_re, R4_im);
    XT_MSUB_S(X4_im, X0_im, R4_re);
    XT_LSXP(R4_im, pR4, -SZ_F32);
    XT_LSXP(R4_re, pR4, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X4_re, X1_re, R4_re);
    XT_MADD_S(X4_re, X1_im, R4_im);
    XT_MSUB_S(X4_im, X1_re, R4_im);
    XT_MSUB_S(X4_im, X1_im, R4_re);
    XT_LSXP(R4_im, pR4, -SZ_F32);
    XT_LSXP(R4_re, pR4, -7 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X4_re, X2_re, R4_re);
    XT_MADD_S(X4_re, X2_im, R4_im);
    XT_MSUB_S(X4_im, X2_re, R4_im);
    XT_MSUB_S(X4_im, X2_im, R4_re);
    XT_LSXP(R4_im, pR4, -SZ_F32);
    XT_LSIP(R4_re, pR4, 0);
    XT_MSUB_S(X4_re, X3_re, R4_re);
    XT_MADD_S(X4_re, X3_im, R4_im);
    XT_MSUB_S(X4_im, X3_re, R4_im);
    XT_MSUB_S(X4_im, X3_im, R4_re);
    X4_re = XT_MUL_S(X4_re, D0_re);
    X4_im = XT_MUL_S(X4_im, D0_im);
    XT_SSXP(X4_im, pX, -SZ_F32);
    XT_SSXP(X4_re, pX, -SZ_F32);

    // X5
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X5_im, pY, -SZ_F32);
    XT_LSXP(X5_re, pY, -SZ_F32);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X5_re, X0_re, R5_re);
    XT_MADD_S(X5_re, X0_im, R5_im);
    XT_MSUB_S(X5_im, X0_re, R5_im);
    XT_MSUB_S(X5_im, X0_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X5_re, X1_re, R5_re);
    XT_MADD_S(X5_re, X1_im, R5_im);
    XT_MSUB_S(X5_im, X1_re, R5_im);
    XT_MSUB_S(X5_im, X1_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -7 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X5_re, X2_re, R5_re);
    XT_MADD_S(X5_re, X2_im, R5_im);
    XT_MSUB_S(X5_im, X2_re, R5_im);
    XT_MSUB_S(X5_im, X2_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -6 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X5_re, X3_re, R5_re);
    XT_MADD_S(X5_re, X3_im, R5_im);
    XT_MSUB_S(X5_im, X3_re, R5_im);
    XT_MSUB_S(X5_im, X3_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSIP(R5_re, pR5, 0);
    XT_MSUB_S(X5_re, X4_re, R5_re);
    XT_MADD_S(X5_re, X4_im, R5_im);
    XT_MSUB_S(X5_im, X4_re, R5_im);
    XT_MSUB_S(X5_im, X4_im, R5_re);
    X5_re = XT_MUL_S(X5_re, D0_re);
    X5_im = XT_MUL_S(X5_im, D0_im);
    XT_SSXP(X5_im, pX, -SZ_F32);
    XT_SSXP(X5_re, pX, -SZ_F32);

    // X6
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X6_im, pY, -SZ_F32);
    XT_LSXP(X6_re, pY, -SZ_F32);
    XT_LSXP(R6_im, pR6, -SZ_F32);
    XT_LSXP(R6_re, pR6, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X6_re, X0_re, R6_re);
    XT_MADD_S(X6_re, X0_im, R6_im);
    XT_MSUB_S(X6_im, X0_re, R6_im);
    XT_MSUB_S(X6_im, X0_im, R6_re);
    XT_LSXP(R6_im, pR6, -SZ_F32);
    XT_LSXP(R6_re, pR6, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X6_re, X1_re, R6_re);
    XT_MADD_S(X6_re, X1_im, R6_im);
    XT_MSUB_S(X6_im, X1_re, R6_im);
    XT_MSUB_S(X6_im, X1_im, R6_re);
    XT_LSXP(R6_im, pR6, -SZ_F32);
    XT_LSXP(R6_re, pR6, -7 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X6_re, X2_re, R6_re);
    XT_MADD_S(X6_re, X2_im, R6_im);
    XT_MSUB_S(X6_im, X2_re, R6_im);
    XT_MSUB_S(X6_im, X2_im, R6_re);
    XT_LSXP(R6_im, pR6, -SZ_F32);
    XT_LSXP(R6_re, pR6, -6 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X6_re, X3_re, R6_re);
    XT_MADD_S(X6_re, X3_im, R6_im);
    XT_MSUB_S(X6_im, X3_re, R6_im);
    XT_MSUB_S(X6_im, X3_im, R6_re);
    XT_LSXP(R6_im, pR6, -SZ_F32);
    XT_LSXP(R6_re, pR6, -5 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X6_re, X4_re, R6_re);
    XT_MADD_S(X6_re, X4_im, R6_im);
    XT_MSUB_S(X6_im, X4_re, R6_im);
    XT_MSUB_S(X6_im, X4_im, R6_re);
    XT_LSXP(R6_im, pR6, -SZ_F32);
    XT_LSIP(R6_re, pR6, 0);
    XT_MSUB_S(X6_re, X5_re, R6_re);
    XT_MADD_S(X6_re, X5_im, R6_im);
    XT_MSUB_S(X6_im, X5_re, R6_im);
    XT_MSUB_S(X6_im, X5_im, R6_re);
    X6_re = XT_MUL_S(X6_re, D0_re);
    X6_im = XT_MUL_S(X6_im, D0_im);
    XT_SSXP(X6_im, pX, -SZ_F32);
    XT_SSXP(X6_re, pX, -SZ_F32);

    // X7
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X7_im, pY, -SZ_F32);
    XT_LSXP(X7_re, pY, -SZ_F32);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSXP(R7_re, pR7, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X7_re, X0_re, R7_re);
    XT_MADD_S(X7_re, X0_im, R7_im);
    XT_MSUB_S(X7_im, X0_re, R7_im);
    XT_MSUB_S(X7_im, X0_im, R7_re);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSXP(R7_re, pR7, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X7_re, X1_re, R7_re);
    XT_MADD_S(X7_re, X1_im, R7_im);
    XT_MSUB_S(X7_im, X1_re, R7_im);
    XT_MSUB_S(X7_im, X1_im, R7_re);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSXP(R7_re, pR7, -7 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X7_re, X2_re, R7_re);
    XT_MADD_S(X7_re, X2_im, R7_im);
    XT_MSUB_S(X7_im, X2_re, R7_im);
    XT_MSUB_S(X7_im, X2_im, R7_re);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSXP(R7_re, pR7, -6 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X7_re, X3_re, R7_re);
    XT_MADD_S(X7_re, X3_im, R7_im);
    XT_MSUB_S(X7_im, X3_re, R7_im);
    XT_MSUB_S(X7_im, X3_im, R7_re);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSXP(R7_re, pR7, -5 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X7_re, X4_re, R7_re);
    XT_MADD_S(X7_re, X4_im, R7_im);
    XT_MSUB_S(X7_im, X4_re, R7_im);
    XT_MSUB_S(X7_im, X4_im, R7_re);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSXP(R7_re, pR7, -4 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X7_re, X5_re, R7_re);
    XT_MADD_S(X7_re, X5_im, R7_im);
    XT_MSUB_S(X7_im, X5_re, R7_im);
    XT_MSUB_S(X7_im, X5_im, R7_re);
    XT_LSXP(R7_im, pR7, -SZ_F32);
    XT_LSIP(R7_re, pR7, 0);
    XT_MSUB_S(X7_re, X6_re, R7_re);
    XT_MADD_S(X7_re, X6_im, R7_im);
    XT_MSUB_S(X7_im, X6_re, R7_im);
    XT_MSUB_S(X7_im, X6_im, R7_re);
    X7_re = XT_MUL_S(X7_re, D0_re);
    X7_im = XT_MUL_S(X7_im, D0_im);
    XT_SSXP(X7_im, pX, -SZ_F32);
    XT_SSXP(X7_re, pX, -SZ_F32);

    // X8
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X8_im, pY, -SZ_F32);
    XT_LSXP(X8_re, pY, -SZ_F32);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X8_re, X0_re, R8_re);
    XT_MADD_S(X8_re, X0_im, R8_im);
    XT_MSUB_S(X8_im, X0_re, R8_im);
    XT_MSUB_S(X8_im, X0_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X8_re, X1_re, R8_re);
    XT_MADD_S(X8_re, X1_im, R8_im);
    XT_MSUB_S(X8_im, X1_re, R8_im);
    XT_MSUB_S(X8_im, X1_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -7 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X8_re, X2_re, R8_re);
    XT_MADD_S(X8_re, X2_im, R8_im);
    XT_MSUB_S(X8_im, X2_re, R8_im);
    XT_MSUB_S(X8_im, X2_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -6 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X8_re, X3_re, R8_re);
    XT_MADD_S(X8_re, X3_im, R8_im);
    XT_MSUB_S(X8_im, X3_re, R8_im);
    XT_MSUB_S(X8_im, X3_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -5 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X8_re, X4_re, R8_re);
    XT_MADD_S(X8_re, X4_im, R8_im);
    XT_MSUB_S(X8_im, X4_re, R8_im);
    XT_MSUB_S(X8_im, X4_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -4 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X8_re, X5_re, R8_re);
    XT_MADD_S(X8_re, X5_im, R8_im);
    XT_MSUB_S(X8_im, X5_re, R8_im);
    XT_MSUB_S(X8_im, X5_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSXP(R8_re, pR8, -3 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X8_re, X6_re, R8_re);
    XT_MADD_S(X8_re, X6_im, R8_im);
    XT_MSUB_S(X8_im, X6_re, R8_im);
    XT_MSUB_S(X8_im, X6_im, R8_re);
    XT_LSXP(R8_im, pR8, -SZ_F32);
    XT_LSIP(R8_re, pR8, 0);
    XT_MSUB_S(X8_re, X7_re, R8_re);
    XT_MADD_S(X8_re, X7_im, R8_im);
    XT_MSUB_S(X8_im, X7_re, R8_im);
    XT_MSUB_S(X8_im, X7_im, R8_re);
    X8_re = XT_MUL_S(X8_re, D0_re);
    X8_im = XT_MUL_S(X8_im, D0_im);
    XT_SSXP(X8_im, pX, -SZ_F32);
    XT_SSXP(X8_re, pX, -SZ_F32);

    // X9
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X9_im, pY, -SZ_F32);
    XT_LSXP(X9_re, pY, -SZ_F32);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -9 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X9_re, X0_re, R9_re);
    XT_MADD_S(X9_re, X0_im, R9_im);
    XT_MSUB_S(X9_im, X0_re, R9_im);
    XT_MSUB_S(X9_im, X0_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -8 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X1_re, R9_re);
    XT_MADD_S(X9_re, X1_im, R9_im);
    XT_MSUB_S(X9_im, X1_re, R9_im);
    XT_MSUB_S(X9_im, X1_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -7 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X2_re, R9_re);
    XT_MADD_S(X9_re, X2_im, R9_im);
    XT_MSUB_S(X9_im, X2_re, R9_im);
    XT_MSUB_S(X9_im, X2_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -6 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X3_re, R9_re);
    XT_MADD_S(X9_re, X3_im, R9_im);
    XT_MSUB_S(X9_im, X3_re, R9_im);
    XT_MSUB_S(X9_im, X3_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -5 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X4_re, R9_re);
    XT_MADD_S(X9_re, X4_im, R9_im);
    XT_MSUB_S(X9_im, X4_re, R9_im);
    XT_MSUB_S(X9_im, X4_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -4 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X5_re, R9_re);
    XT_MADD_S(X9_re, X5_im, R9_im);
    XT_MSUB_S(X9_im, X5_re, R9_im);
    XT_MSUB_S(X9_im, X5_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -3 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X6_re, R9_re);
    XT_MADD_S(X9_re, X6_im, R9_im);
    XT_MSUB_S(X9_im, X6_re, R9_im);
    XT_MSUB_S(X9_im, X6_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSXP(R9_re, pR9, -2 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X9_re, X7_re, R9_re);
    XT_MADD_S(X9_re, X7_im, R9_im);
    XT_MSUB_S(X9_im, X7_re, R9_im);
    XT_MSUB_S(X9_im, X7_im, R9_re);
    XT_LSXP(R9_im, pR9, -SZ_F32);
    XT_LSIP(R9_re, pR9, 0);
    XT_MSUB_S(X9_re, X8_re, R9_re);
    XT_MADD_S(X9_re, X8_im, R9_im);
    XT_MSUB_S(X9_im, X8_re, R9_im);
    XT_MSUB_S(X9_im, X8_im, R9_re);
    X9_re = XT_MUL_S(X9_re, D0_re);
    X9_im = XT_MUL_S(X9_im, D0_im);
    XT_SSXP(X9_im, pX, -SZ_F32);
    XT_SSXP(X9_re, pX, -SZ_F32);
}

size_t  cmatcholbkwsubst10x10f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, cmatcholbkwsubst10x10f, (void * pScr,
    complex_float * x,
    const complex_float * R,
    const complex_float * D,
    const complex_float * y))

size_t  cmatcholbkwsubst10x10f_getScratchSize()
{
    return 0;
}
#endif
