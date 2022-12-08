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
#include "common_fpu.h"
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*--------------------------------------------------------
Forward recursion P=1
Input:
R[((N+1)*N)/2]  upper triangular matrix R
stride          width of matrix A'*B
Z[N*1]          column in matrix A'*B
D[N]            reciprocals of main diagonal
Output:
y[N*1]		Decision matrix y
N = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
--------------------------------------------------------*/
void cplxcholFwdrec8f(complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * Z, int stride)
{
    xtfloatx2 * restrict pY = (xtfloatx2 *)y;
    const xtfloatx2 * restrict pZ = (const xtfloatx2 *)Z;
    const xtfloatx2 * restrict pR1 = (const xtfloatx2 *)(R + 1);/* + ((n*(n+1))/2) */
    const xtfloatx2 * restrict pR2 = (const xtfloatx2 *)(R + 3);
    const xtfloatx2 * restrict pR3 = (const xtfloatx2 *)(R + 6);
    const xtfloatx2 * restrict pR4 = (const xtfloatx2 *)(R + 10);
    const xtfloatx2 * restrict pR5 = (const xtfloatx2 *)(R + 15);
    const xtfloatx2 * restrict pR6 = (const xtfloatx2 *)(R + 21);
    const xtfloatx2 * restrict pR7 = (const xtfloatx2 *)(R + 28);
    const xtfloatx2 * restrict pD = (const xtfloatx2 *)D;
    xtfloatx2 Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7;
    xtfloatx2 Z1t, Z2t, Z3t, Z4t, Z5t, Z6t, Z7t;
    xtfloatx2 D0, D1, D2, D3, D4, D5, D6, D7;
    xtfloatx2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    xtfloatx2 R0, R1, R2, R3, R4, R5, R6;

    Z1t = Z2t = Z3t = Z4t = Z5t = Z6t = Z7t = XT_CONST_S(0);
    // Y0
    XT_LSX2XP(Z0, pZ, stride*SZ_CF32);
    XT_LSX2IP(D0, pD, SZ_CF32);
    Y0 = XT_MUL_SX2(Z0, D0);
    XT_SSX2IP(Y0, pY, SZ_CF32);

    // Y1
    XT_LSX2XP(Z1, pZ, stride*SZ_CF32);
    XT_LSX2IP(D1, pD, SZ_CF32);
    XT_LSX2IP(R0, pR1, 0);
    XT_MADDMUX_S(Z1, R0, Y0, 2);
    XT_MADDMUX_S(Z1t, R0, Y0, 1);
    Z1 = XT_ADD_SX2(Z1, Z1t);
    Y1 = XT_MUL_SX2(Z1, D1);
    XT_SSX2IP(Y1, pY, SZ_CF32);

    // Y2
    XT_LSX2XP(Z2, pZ, stride*SZ_CF32);
    XT_LSX2IP(D2, pD, SZ_CF32);
    XT_LSX2IP(R0, pR2, SZ_CF32);
    XT_LSX2IP(R1, pR2, 0);
    XT_MADDMUX_S(Z2, R0, Y0, 2);
    XT_MADDMUX_S(Z2t, R0, Y0, 1);
    XT_MADDMUX_S(Z2, R1, Y1, 2);
    XT_MADDMUX_S(Z2t, R1, Y1, 1);
    Z2 = XT_ADD_SX2(Z2, Z2t);
    Y2 = XT_MUL_SX2(Z2, D2);
    XT_SSX2IP(Y2, pY, SZ_CF32);

    // Y3
    XT_LSX2XP(Z3, pZ, stride*SZ_CF32);
    XT_LSX2IP(D3, pD, SZ_CF32);
    XT_LSX2IP(R0, pR3, SZ_CF32);
    XT_LSX2IP(R1, pR3, SZ_CF32);
    XT_LSX2IP(R2, pR3, 0);
    XT_MADDMUX_S(Z3, R0, Y0, 2);
    XT_MADDMUX_S(Z3t, R0, Y0, 1);
    XT_MADDMUX_S(Z3, R1, Y1, 2);
    XT_MADDMUX_S(Z3t, R1, Y1, 1);
    XT_MADDMUX_S(Z3, R2, Y2, 2);
    XT_MADDMUX_S(Z3t, R2, Y2, 1);
    Z3 = XT_ADD_SX2(Z3, Z3t);
    Y3 = XT_MUL_SX2(Z3, D3);
    XT_SSX2IP(Y3, pY, SZ_CF32);

    // Y4
    XT_LSX2XP(Z4, pZ, stride*SZ_CF32);
    XT_LSX2IP(D4, pD, SZ_CF32);
    XT_LSX2IP(R0, pR4, SZ_CF32);
    XT_LSX2IP(R1, pR4, SZ_CF32);
    XT_LSX2IP(R2, pR4, SZ_CF32);
    XT_LSX2IP(R3, pR4, 0);
    XT_MADDMUX_S(Z4, R0, Y0, 2);
    XT_MADDMUX_S(Z4t, R0, Y0, 1);
    XT_MADDMUX_S(Z4, R1, Y1, 2);
    XT_MADDMUX_S(Z4t, R1, Y1, 1);
    XT_MADDMUX_S(Z4, R2, Y2, 2);
    XT_MADDMUX_S(Z4t, R2, Y2, 1);
    XT_MADDMUX_S(Z4, R3, Y3, 2);
    XT_MADDMUX_S(Z4t, R3, Y3, 1);
    Z4 = XT_ADD_SX2(Z4, Z4t);
    Y4 = XT_MUL_SX2(Z4, D4);
    XT_SSX2IP(Y4, pY, SZ_CF32);

    // Y5
    XT_LSX2XP(Z5, pZ, stride*SZ_CF32);
    XT_LSX2IP(D5, pD, SZ_CF32);
    XT_LSX2IP(R0, pR5, SZ_CF32);
    XT_LSX2IP(R1, pR5, SZ_CF32);
    XT_LSX2IP(R2, pR5, SZ_CF32);
    XT_LSX2IP(R3, pR5, SZ_CF32);
    XT_LSX2IP(R4, pR5, 0);
    XT_MADDMUX_S(Z5, R0, Y0, 2);
    XT_MADDMUX_S(Z5t, R0, Y0, 1);
    XT_MADDMUX_S(Z5, R1, Y1, 2);
    XT_MADDMUX_S(Z5t, R1, Y1, 1);
    XT_MADDMUX_S(Z5, R2, Y2, 2);
    XT_MADDMUX_S(Z5t, R2, Y2, 1);
    XT_MADDMUX_S(Z5, R3, Y3, 2);
    XT_MADDMUX_S(Z5t, R3, Y3, 1);
    XT_MADDMUX_S(Z5, R4, Y4, 2);
    XT_MADDMUX_S(Z5t, R4, Y4, 1);
    Z5 = XT_ADD_SX2(Z5, Z5t);
    Y5 = XT_MUL_SX2(Z5, D5);
    XT_SSX2IP(Y5, pY, SZ_CF32);

    // Y6
    XT_LSX2XP(Z6, pZ, stride*SZ_CF32);
    XT_LSX2IP(D6, pD, SZ_CF32);
    XT_LSX2IP(R0, pR6, SZ_CF32);
    XT_LSX2IP(R1, pR6, SZ_CF32);
    XT_LSX2IP(R2, pR6, SZ_CF32);
    XT_LSX2IP(R3, pR6, SZ_CF32);
    XT_LSX2IP(R4, pR6, SZ_CF32);
    XT_LSX2IP(R5, pR6, 0);
    XT_MADDMUX_S(Z6, R0, Y0, 2);
    XT_MADDMUX_S(Z6t, R0, Y0, 1);
    XT_MADDMUX_S(Z6, R1, Y1, 2);
    XT_MADDMUX_S(Z6t, R1, Y1, 1);
    XT_MADDMUX_S(Z6, R2, Y2, 2);
    XT_MADDMUX_S(Z6t, R2, Y2, 1);
    XT_MADDMUX_S(Z6, R3, Y3, 2);
    XT_MADDMUX_S(Z6t, R3, Y3, 1);
    XT_MADDMUX_S(Z6, R4, Y4, 2);
    XT_MADDMUX_S(Z6t, R4, Y4, 1);
    XT_MADDMUX_S(Z6, R5, Y5, 2);
    XT_MADDMUX_S(Z6t, R5, Y5, 1);
    Z6 = XT_ADD_SX2(Z6, Z6t);
    Y6 = XT_MUL_SX2(Z6, D6);
    XT_SSX2IP(Y6, pY, SZ_CF32);

    // Y7
    XT_LSX2XP(Z7, pZ, stride*SZ_CF32);
    XT_LSX2IP(D7, pD, SZ_CF32);
    XT_LSX2IP(R0, pR7, SZ_CF32);
    XT_LSX2IP(R1, pR7, SZ_CF32);
    XT_LSX2IP(R2, pR7, SZ_CF32);
    XT_LSX2IP(R3, pR7, SZ_CF32);
    XT_LSX2IP(R4, pR7, SZ_CF32);
    XT_LSX2IP(R5, pR7, SZ_CF32);
    XT_LSX2IP(R6, pR7, 0);
    XT_MADDMUX_S(Z7, R0, Y0, 2);
    XT_MADDMUX_S(Z7t, R0, Y0, 1);
    XT_MADDMUX_S(Z7, R1, Y1, 2);
    XT_MADDMUX_S(Z7t, R1, Y1, 1);
    XT_MADDMUX_S(Z7, R2, Y2, 2);
    XT_MADDMUX_S(Z7t, R2, Y2, 1);
    XT_MADDMUX_S(Z7, R3, Y3, 2);
    XT_MADDMUX_S(Z7t, R3, Y3, 1);
    XT_MADDMUX_S(Z7, R4, Y4, 2);
    XT_MADDMUX_S(Z7t, R4, Y4, 1);
    XT_MADDMUX_S(Z7, R5, Y5, 2);
    XT_MADDMUX_S(Z7t, R5, Y5, 1);
    XT_MADDMUX_S(Z7, R6, Y6, 2);
    XT_MADDMUX_S(Z7t, R6, Y6, 1);
    Z7 = XT_ADD_SX2(Z7, Z7t);
    Y7 = XT_MUL_SX2(Z7, D7);
    XT_SSX2IP(Y7, pY, SZ_CF32);
}
#elif (HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*--------------------------------------------------------
Forward recursion P=1
Input:
R[((N+1)*N)/2]  upper triangular matrix R
stride          width of matrix A'*B
Z[N*1]          column in matrix A'*B
D[N]            reciprocals of main diagonal
Output:
y[N*1]		Decision matrix y
N = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
--------------------------------------------------------*/
void cplxcholFwdrec8f(complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * Z, int stride)
{
    xtfloat * restrict pY = (xtfloat*)y;
    const xtfloat * restrict pZ = (const xtfloat*)Z;
    const xtfloat * restrict pD = (const xtfloat*)D;
    const xtfloat * restrict pR = (const xtfloat*)(R + 1);
    xtfloat Y0_re, Y1_re, Y2_re, Y3_re, Y4_re, Y5_re, Y6_re, Y7_re;
    xtfloat Y0_im, Y1_im, Y2_im, Y3_im, Y4_im, Y5_im, Y6_im, Y7_im;
    xtfloat R_re, R_im;
    xtfloat D_re, D_im;

    //Y0
    XT_LSIP(Y0_re, pZ, SZ_F32);
    XT_LSXP(Y0_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    Y0_re = XT_MUL_S(Y0_re, D_re);
    Y0_im = XT_MUL_S(Y0_im, D_im);
    XT_SSIP(Y0_re, pY, SZ_F32);
    XT_SSIP(Y0_im, pY, SZ_F32);

    //Y1
    XT_LSIP(Y1_re, pZ, SZ_F32);
    XT_LSXP(Y1_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y1_re, Y0_re, R_re);
    XT_MSUB_S(Y1_re, Y0_im, R_im);
    XT_MSUB_S(Y1_im, Y0_im, R_re);
    XT_MADD_S(Y1_im, Y0_re, R_im);

    Y1_re = XT_MUL_S(Y1_re, D_re);
    Y1_im = XT_MUL_S(Y1_im, D_im);
    XT_SSIP(Y1_re, pY, SZ_F32);
    XT_SSIP(Y1_im, pY, SZ_F32);

    //Y2
    XT_LSIP(Y2_re, pZ, SZ_F32);
    XT_LSXP(Y2_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y2_re, Y0_re, R_re);
    XT_MSUB_S(Y2_re, Y0_im, R_im);
    XT_MSUB_S(Y2_im, Y0_im, R_re);
    XT_MADD_S(Y2_im, Y0_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y2_re, Y1_re, R_re);
    XT_MSUB_S(Y2_re, Y1_im, R_im);
    XT_MSUB_S(Y2_im, Y1_im, R_re);
    XT_MADD_S(Y2_im, Y1_re, R_im);

    Y2_re = XT_MUL_S(Y2_re, D_re);
    Y2_im = XT_MUL_S(Y2_im, D_im);
    XT_SSIP(Y2_re, pY, SZ_F32);
    XT_SSIP(Y2_im, pY, SZ_F32);

    //Y3
    XT_LSIP(Y3_re, pZ, SZ_F32);
    XT_LSXP(Y3_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y3_re, Y0_re, R_re);
    XT_MSUB_S(Y3_re, Y0_im, R_im);
    XT_MSUB_S(Y3_im, Y0_im, R_re);
    XT_MADD_S(Y3_im, Y0_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y3_re, Y1_re, R_re);
    XT_MSUB_S(Y3_re, Y1_im, R_im);
    XT_MSUB_S(Y3_im, Y1_im, R_re);
    XT_MADD_S(Y3_im, Y1_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y3_re, Y2_re, R_re);
    XT_MSUB_S(Y3_re, Y2_im, R_im);
    XT_MSUB_S(Y3_im, Y2_im, R_re);
    XT_MADD_S(Y3_im, Y2_re, R_im);
    Y3_re = XT_MUL_S(Y3_re, D_re);
    Y3_im = XT_MUL_S(Y3_im, D_im);

    XT_SSIP(Y3_re, pY, SZ_F32);
    XT_SSIP(Y3_im, pY, SZ_F32);

    //Y4
    XT_LSIP(Y4_re, pZ, SZ_F32);
    XT_LSXP(Y4_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y4_re, Y0_re, R_re);
    XT_MSUB_S(Y4_re, Y0_im, R_im);
    XT_MSUB_S(Y4_im, Y0_im, R_re);
    XT_MADD_S(Y4_im, Y0_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y4_re, Y1_re, R_re);
    XT_MSUB_S(Y4_re, Y1_im, R_im);
    XT_MSUB_S(Y4_im, Y1_im, R_re);
    XT_MADD_S(Y4_im, Y1_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y4_re, Y2_re, R_re);
    XT_MSUB_S(Y4_re, Y2_im, R_im);
    XT_MSUB_S(Y4_im, Y2_im, R_re);
    XT_MADD_S(Y4_im, Y2_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y4_re, Y3_re, R_re);
    XT_MSUB_S(Y4_re, Y3_im, R_im);
    XT_MSUB_S(Y4_im, Y3_im, R_re);
    XT_MADD_S(Y4_im, Y3_re, R_im);
    Y4_re = XT_MUL_S(Y4_re, D_re);
    Y4_im = XT_MUL_S(Y4_im, D_im);

    XT_SSIP(Y4_re, pY, SZ_F32);
    XT_SSIP(Y4_im, pY, SZ_F32);

    //Y5
    XT_LSIP(Y5_re, pZ, SZ_F32);
    XT_LSXP(Y5_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y5_re, Y0_re, R_re);
    XT_MSUB_S(Y5_re, Y0_im, R_im);
    XT_MSUB_S(Y5_im, Y0_im, R_re);
    XT_MADD_S(Y5_im, Y0_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y5_re, Y1_re, R_re);
    XT_MSUB_S(Y5_re, Y1_im, R_im);
    XT_MSUB_S(Y5_im, Y1_im, R_re);
    XT_MADD_S(Y5_im, Y1_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y5_re, Y2_re, R_re);
    XT_MSUB_S(Y5_re, Y2_im, R_im);
    XT_MSUB_S(Y5_im, Y2_im, R_re);
    XT_MADD_S(Y5_im, Y2_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y5_re, Y3_re, R_re);
    XT_MSUB_S(Y5_re, Y3_im, R_im);
    XT_MSUB_S(Y5_im, Y3_im, R_re);
    XT_MADD_S(Y5_im, Y3_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y5_re, Y4_re, R_re);
    XT_MSUB_S(Y5_re, Y4_im, R_im);
    XT_MSUB_S(Y5_im, Y4_im, R_re);
    XT_MADD_S(Y5_im, Y4_re, R_im);
    Y5_re = XT_MUL_S(Y5_re, D_re);
    Y5_im = XT_MUL_S(Y5_im, D_im);

    XT_SSIP(Y5_re, pY, SZ_F32);
    XT_SSIP(Y5_im, pY, SZ_F32);

    //Y6
    XT_LSIP(Y6_re, pZ, SZ_F32);
    XT_LSXP(Y6_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y6_re, Y0_re, R_re);
    XT_MSUB_S(Y6_re, Y0_im, R_im);
    XT_MSUB_S(Y6_im, Y0_im, R_re);
    XT_MADD_S(Y6_im, Y0_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y6_re, Y1_re, R_re);
    XT_MSUB_S(Y6_re, Y1_im, R_im);
    XT_MSUB_S(Y6_im, Y1_im, R_re);
    XT_MADD_S(Y6_im, Y1_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y6_re, Y2_re, R_re);
    XT_MSUB_S(Y6_re, Y2_im, R_im);
    XT_MSUB_S(Y6_im, Y2_im, R_re);
    XT_MADD_S(Y6_im, Y2_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y6_re, Y3_re, R_re);
    XT_MSUB_S(Y6_re, Y3_im, R_im);
    XT_MSUB_S(Y6_im, Y3_im, R_re);
    XT_MADD_S(Y6_im, Y3_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y6_re, Y4_re, R_re);
    XT_MSUB_S(Y6_re, Y4_im, R_im);
    XT_MSUB_S(Y6_im, Y4_im, R_re);
    XT_MADD_S(Y6_im, Y4_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y6_re, Y5_re, R_re);
    XT_MSUB_S(Y6_re, Y5_im, R_im);
    XT_MSUB_S(Y6_im, Y5_im, R_re);
    XT_MADD_S(Y6_im, Y5_re, R_im);
    Y6_re = XT_MUL_S(Y6_re, D_re);
    Y6_im = XT_MUL_S(Y6_im, D_im);

    XT_SSIP(Y6_re, pY, SZ_F32);
    XT_SSIP(Y6_im, pY, SZ_F32);

    //Y7
    XT_LSIP(Y7_re, pZ, SZ_F32);
    XT_LSXP(Y7_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y7_re, Y0_re, R_re);
    XT_MSUB_S(Y7_re, Y0_im, R_im);
    XT_MSUB_S(Y7_im, Y0_im, R_re);
    XT_MADD_S(Y7_im, Y0_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y7_re, Y1_re, R_re);
    XT_MSUB_S(Y7_re, Y1_im, R_im);
    XT_MSUB_S(Y7_im, Y1_im, R_re);
    XT_MADD_S(Y7_im, Y1_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y7_re, Y2_re, R_re);
    XT_MSUB_S(Y7_re, Y2_im, R_im);
    XT_MSUB_S(Y7_im, Y2_im, R_re);
    XT_MADD_S(Y7_im, Y2_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y7_re, Y3_re, R_re);
    XT_MSUB_S(Y7_re, Y3_im, R_im);
    XT_MSUB_S(Y7_im, Y3_im, R_re);
    XT_MADD_S(Y7_im, Y3_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y7_re, Y4_re, R_re);
    XT_MSUB_S(Y7_re, Y4_im, R_im);
    XT_MSUB_S(Y7_im, Y4_im, R_re);
    XT_MADD_S(Y7_im, Y4_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, SZ_F32);
    XT_MSUB_S(Y7_re, Y5_re, R_re);
    XT_MSUB_S(Y7_re, Y5_im, R_im);
    XT_MSUB_S(Y7_im, Y5_im, R_re);
    XT_MADD_S(Y7_im, Y5_re, R_im);
    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y7_re, Y6_re, R_re);
    XT_MSUB_S(Y7_re, Y6_im, R_im);
    XT_MSUB_S(Y7_im, Y6_im, R_re);
    XT_MADD_S(Y7_im, Y6_re, R_im);
    Y7_re = XT_MUL_S(Y7_re, D_re);
    Y7_im = XT_MUL_S(Y7_im, D_im);

    XT_SSIP(Y7_re, pY, SZ_F32);
    XT_SSIP(Y7_im, pY, SZ_F32);
}
#endif
