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
void  cmatcholbkwsubst6x6f(void * pScr,
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

    const xtfloatx2 * restrict pR5 = (const xtfloatx2*)(R + (6 * 5) / 2 + 0); //N*(N-1)/2
    const xtfloatx2 * restrict pR4 = (const xtfloatx2*)(R + (6 * 5) / 2 + 1);
    const xtfloatx2 * restrict pR3 = (const xtfloatx2*)(R + (6 * 5) / 2 + 2);
    const xtfloatx2 * restrict pR2 = (const xtfloatx2*)(R + (6 * 5) / 2 + 3);
    const xtfloatx2 * restrict pR1 = (const xtfloatx2*)(R + (6 * 5) / 2 + 4);
    const xtfloatx2 * restrict pY = (const xtfloatx2*)(y + 5);
    const xtfloatx2 * restrict pD = (const xtfloatx2*)(D + 5);
    xtfloatx2 * restrict pX = (xtfloatx2*)(x + 5);
    xtfloatx2 D0;
    xtfloatx2 R1, R2, R3, R4, R5;
    xtfloatx2 X0, X1, X2, X3, X4, X5, Xt;

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
    XT_LSX2XP(R2, pR2, -SZ_CF32 * 5); // N-1
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
    XT_LSX2XP(R3, pR3, -SZ_CF32 * 5);
    XT_MADDMUX_S(X3, X0, R3, 2);
    XT_MADDMUX_S(Xt, X0, R3, 3);
    XT_LSX2XP(R3, pR3, -SZ_CF32 * 4);
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
    XT_LSX2XP(R4, pR4, -SZ_CF32 * 5);
    XT_MADDMUX_S(X4, X0, R4, 2);
    XT_MADDMUX_S(Xt, X0, R4, 3);
    XT_LSX2XP(R4, pR4, -SZ_CF32 * 4);
    XT_MADDMUX_S(X4, X1, R4, 2);
    XT_MADDMUX_S(Xt, X1, R4, 3);
    XT_LSX2XP(R4, pR4, -SZ_CF32 * 3);
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
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 5);
    XT_MADDMUX_S(X5, X0, R5, 2);
    XT_MADDMUX_S(Xt, X0, R5, 3);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 4);
    XT_MADDMUX_S(X5, X1, R5, 2);
    XT_MADDMUX_S(Xt, X1, R5, 3);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 3);
    XT_MADDMUX_S(X5, X2, R5, 2);
    XT_MADDMUX_S(Xt, X2, R5, 3);
    XT_LSX2XP(R5, pR5, -SZ_CF32 * 2);
    XT_MADDMUX_S(X5, X3, R5, 2);
    XT_MADDMUX_S(Xt, X3, R5, 3);
    XT_LSX2IP(R5, pR5, 0);
    XT_MADDMUX_S(X5, X4, R5, 2);
    XT_MADDMUX_S(Xt, X4, R5, 3);
    X5 = XT_ADD_SX2(X5, Xt);
    X5 = XT_MUL_SX2(X5, D0);
    XT_SSX2XP(X5, pX, -SZ_CF32);

}

size_t  cmatcholbkwsubst6x6f_getScratchSize()
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
void  cmatcholbkwsubst6x6f(void * pScr,
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

    const xtfloat * restrict pR5 = (const xtfloat*)(R + (6 * 5) / 2 + 0) + 1; //N*(N-1)/2
    const xtfloat * restrict pR4 = (const xtfloat*)(R + (6 * 5) / 2 + 1) + 1;
    const xtfloat * restrict pR3 = (const xtfloat*)(R + (6 * 5) / 2 + 2) + 1;
    const xtfloat * restrict pR2 = (const xtfloat*)(R + (6 * 5) / 2 + 3) + 1;
    const xtfloat * restrict pR1 = (const xtfloat*)(R + (6 * 5) / 2 + 4) + 1;
    const xtfloat * restrict pY = (const xtfloat*)(y + 5) + 1;
    const xtfloat * restrict pD = (const xtfloat*)(D + 5) + 1;
    xtfloat * restrict pX = (xtfloat*)(x + 5) + 1;
    xtfloat D0_re, D0_im;
    xtfloat R1_re, R2_re, R3_re, R4_re, R5_re;
    xtfloat R1_im, R2_im, R3_im, R4_im, R5_im;
    xtfloat X0_re, X1_re, X2_re, X3_re, X4_re, X5_re;
    xtfloat X0_im, X1_im, X2_im, X3_im, X4_im, X5_im;

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
    XT_LSXP(R2_re, pR2, -5 * SZ_CF32 + SZ_F32); //N-1
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
    XT_LSXP(R3_re, pR3, -5 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X3_re, X0_re, R3_re);
    XT_MADD_S(X3_re, X0_im, R3_im);
    XT_MSUB_S(X3_im, X0_re, R3_im);
    XT_MSUB_S(X3_im, X0_im, R3_re);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSXP(R3_re, pR3, -4 * SZ_CF32 + SZ_F32);
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
    XT_LSXP(R4_re, pR4, -5 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X4_re, X0_re, R4_re);
    XT_MADD_S(X4_re, X0_im, R4_im);
    XT_MSUB_S(X4_im, X0_re, R4_im);
    XT_MSUB_S(X4_im, X0_im, R4_re);
    XT_LSXP(R4_im, pR4, -SZ_F32);
    XT_LSXP(R4_re, pR4, -4 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X4_re, X1_re, R4_re);
    XT_MADD_S(X4_re, X1_im, R4_im);
    XT_MSUB_S(X4_im, X1_re, R4_im);
    XT_MSUB_S(X4_im, X1_im, R4_re);
    XT_LSXP(R4_im, pR4, -SZ_F32 );
    XT_LSXP(R4_re, pR4, -3 * SZ_CF32 + SZ_F32);
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
    XT_LSXP(R5_re, pR5, -5 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X5_re, X0_re, R5_re);
    XT_MADD_S(X5_re, X0_im, R5_im);
    XT_MSUB_S(X5_im, X0_re, R5_im);
    XT_MSUB_S(X5_im, X0_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -4 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X5_re, X1_re, R5_re);
    XT_MADD_S(X5_re, X1_im, R5_im);
    XT_MSUB_S(X5_im, X1_re, R5_im);
    XT_MSUB_S(X5_im, X1_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -3 * SZ_CF32 + SZ_F32);
    XT_MSUB_S(X5_re, X2_re, R5_re);
    XT_MADD_S(X5_re, X2_im, R5_im);
    XT_MSUB_S(X5_im, X2_re, R5_im);
    XT_MSUB_S(X5_im, X2_im, R5_re);
    XT_LSXP(R5_im, pR5, -SZ_F32);
    XT_LSXP(R5_re, pR5, -2 * SZ_CF32 + SZ_F32);
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
}

size_t  cmatcholbkwsubst6x6f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, cmatcholbkwsubst6x6f, (void * pScr,
    complex_float * x,
    const complex_float * R,
    const complex_float * D,
    const complex_float * y))

size_t  cmatcholbkwsubst6x6f_getScratchSize()
{
    return 0;
}
#endif
