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
#define F32_WIDTH (HIFI_SIMD_WIDTH/SZ_F32)
/*
Forward recursion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A'*B
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
Temporary:
pScr        Scratch memory
*/
static void rfwd10x10f(
    float32_t* restrict y,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict Z)
{
    int n, m;
    xtfloat R0, D0;
    xtfloatx2 Yr01, Yr23, Yr45, Yr67, Yr89;
    xtfloat Yw0, Yw1, Yw2, Yw3, Yw4, Yw5, Yw6, Yw7, Yw8, Yw9;
    xtfloatx2 Yw01, Yw23, Yw45, Yw67, Yw89;
    xtfloatx2 * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)R;
    const xtfloat * restrict pD = (const xtfloat *)D;
    const xtfloat * restrict pA = (const xtfloat *)Z;
    for (n = 0; n<10; n++)
    {
        pY = (xtfloatx2*)y;
        XT_LSXP(Yw0, pA, 10 * SZ_F32);
        XT_LSXP(Yw1, pA, 10 * SZ_F32);
        XT_LSXP(Yw2, pA, 10 * SZ_F32);
        XT_LSXP(Yw3, pA, 10 * SZ_F32);
        XT_LSXP(Yw4, pA, 10 * SZ_F32);
        XT_LSXP(Yw5, pA, 10 * SZ_F32);
        XT_LSXP(Yw6, pA, 10 * SZ_F32);
        XT_LSXP(Yw7, pA, 10 * SZ_F32);
        XT_LSXP(Yw8, pA, 10 * SZ_F32);
        XT_LSXP(Yw9, pA, (-10 * (10 - 1) + 1)*SZ_F32);
        Yw01 = XT_SEL32_HH_SX2(Yw0, Yw1);
        Yw23 = XT_SEL32_HH_SX2(Yw2, Yw3);
        Yw45 = XT_SEL32_HH_SX2(Yw4, Yw5);
        Yw67 = XT_SEL32_HH_SX2(Yw6, Yw7);
        Yw89 = XT_SEL32_HH_SX2(Yw8, Yw9);

        for (m = 0; m<n; m++)
        {
            XT_LSX2IP(Yr01, pY, 2 * SZ_F32);
            XT_LSX2IP(Yr23, pY, 2 * SZ_F32);
            XT_LSX2IP(Yr45, pY, 2 * SZ_F32);
            XT_LSX2IP(Yr67, pY, 2 * SZ_F32);
            XT_LSX2IP(Yr89, pY, 2 * SZ_F32);
            XT_LSIP(R0, pR, SZ_F32);

            XT_MSUB_SX2(Yw01, R0, Yr01);
            XT_MSUB_SX2(Yw23, R0, Yr23);
            XT_MSUB_SX2(Yw45, R0, Yr45);
            XT_MSUB_SX2(Yw67, R0, Yr67);
            XT_MSUB_SX2(Yw89, R0, Yr89);
        }

        XT_LSIP(D0, pD, SZ_F32);
        Yw01 = XT_MUL_SX2(Yw01, D0);
        Yw23 = XT_MUL_SX2(Yw23, D0);
        Yw45 = XT_MUL_SX2(Yw45, D0);
        Yw67 = XT_MUL_SX2(Yw67, D0);
        Yw89 = XT_MUL_SX2(Yw89, D0);

        XT_SSX2IP(Yw01, pY, 2 * SZ_F32);
        XT_SSX2IP(Yw23, pY, 2 * SZ_F32);
        XT_SSX2IP(Yw45, pY, 2 * SZ_F32);
        XT_SSX2IP(Yw67, pY, 2 * SZ_F32);
        XT_SSX2IP(Yw89, pY, 2 * SZ_F32);
        pR += 1;
    }
}
/*
backward recursion
does not require tranformed R, reverse inner loop
*/
static void rbkw10x10f(
    float32_t* restrict x,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k;
    xtfloatx2 * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 10 * (10 + 1) / 2 - 1); //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (10 - 1));
    const xtfloatx2 * restrict pY = (xtfloatx2*)(y + (10 * 10 - 2));
    xtfloatx2 Xw01, Xw23, Xw45, Xw67, Xw89;
    xtfloatx2 Xr01, Xr23, Xr45, Xr67, Xr89;
    xtfloat R0, D0;

    for (k = 10 - 1; k >= 0; k--)
    {
        {
            pX = (xtfloatx2*)(x + 10 * 10 - 2); // last 2 elements in X
            pR = (xtfloat*)(tpR--); //points to the end of row k in R
            XT_LSX2XP(Xw89, pY, -2 * SZ_F32);
            XT_LSX2XP(Xw67, pY, -2 * SZ_F32);
            XT_LSX2XP(Xw45, pY, -2 * SZ_F32);
            XT_LSX2XP(Xw23, pY, -2 * SZ_F32);
            XT_LSX2XP(Xw01, pY, -2 * SZ_F32);

            for (m = 0; m < 10 - k - 1; m++)
            {
                XT_LSX2XP(Xr89, pX, -2 * SZ_F32);
                XT_LSX2XP(Xr67, pX, -2 * SZ_F32);
                XT_LSX2XP(Xr45, pX, -2 * SZ_F32);
                XT_LSX2XP(Xr23, pX, -2 * SZ_F32);
                XT_LSX2XP(Xr01, pX, -2 * SZ_F32);

                XT_LSXP(R0, pR, -(10 - 1 - m)*SZ_F32);
                XT_MSUB_SX2(Xw89, Xr89, R0);
                XT_MSUB_SX2(Xw67, Xr67, R0);
                XT_MSUB_SX2(Xw45, Xr45, R0);
                XT_MSUB_SX2(Xw23, Xr23, R0);
                XT_MSUB_SX2(Xw01, Xr01, R0);
            }
            XT_LSXP(D0, pD, -SZ_F32);
            Xw89 = XT_MUL_SX2(Xw89, D0);
            Xw67 = XT_MUL_SX2(Xw67, D0);
            Xw45 = XT_MUL_SX2(Xw45, D0);
            Xw23 = XT_MUL_SX2(Xw23, D0);
            Xw01 = XT_MUL_SX2(Xw01, D0);
            XT_SSX2XP(Xw89, pX, -2 * SZ_F32);
            XT_SSX2XP(Xw67, pX, -2 * SZ_F32);
            XT_SSX2XP(Xw45, pX, -2 * SZ_F32);
            XT_SSX2XP(Xw23, pX, -2 * SZ_F32);
            XT_SSX2XP(Xw01, pX, -2 * SZ_F32);

        }
    }
}
/*-------------------------------------------------------------------------
Matrix (Pseudo) Inversion
Obtain Left Inverse of a matrix using Cholesky Decomposition
The result is matrix x = A^-1
Fixed point API requires explicit setting of fixed point representation of 
input/output matrices as well as for internal temporary matrices such as R 
(Cholesky decomposition) and y (decision of R'y=(A'*B))


Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]      matrix A, for fixed point API, the representation is Q(qA)
sigma2      regularization term, for fixed point API, the 
            representation is Q(2*qA-30)
qRA         qR-qA; difference between fixed point representations of R
            and A matrices (for the fixed point API only). Should be 
            equal or less than 0 (typically -2).
qYBRA       qY-qB+qR-qA, combination of fixed point representations of 
            matrices y, B, R and A (for the fixed point API only). Since 
            for matrix inversion we simply use identity matrix B, we may 
            always suppose qB=31 
qXYR        combination of fixed point representation (matrices R, x and y) 
            qXYR=qX-qY+qR (for the fixed point API only)
Output:
x[N*M]      Left Inverse of the matrix A, for fixed point API, the 
            representation is Q(qX)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  matcholpseudoinv10x10f(void* pScr,
    float32_t *x,
    const float32_t * A,
    const float32_t  sigma2)
{
    float32_t * restrict D; int SD;
    float32_t * restrict R; int SR;
    float32_t * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);
    R = (float32_t *)pScr;
    D = (float32_t *)(((uintptr_t)R) + SR);
    y = (float32_t *)(((uintptr_t)D) + SD);
    pScr = (float32_t *)(((uintptr_t)y) + SY);
    matcholdecomp10x10f(pScr, R, D, A, sigma2);
    rfwd10x10f(y, R, D, A);
    rbkw10x10f(x, R, D, y);
}

size_t   matcholpseudoinv10x10f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);

    s_dc = matcholdecomp10x10f_getScratchSize();

    return SD + SR + SY + s_dc;
}
#elif (HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define F32_WIDTH (HIFI_SIMD_WIDTH/SZ_F32)
/*
Forward recursion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A'*B
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
Temporary:
pScr        Scratch memory
*/
static void rfwd10x10f(
    float32_t* restrict y,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict Z)
{
#if 0
    int n, m, p;
    xtfloat Yw, Yr, R0, D0;
    xtfloat * restrict pY;
    const xtfloat * restrict pR;
    const xtfloat * restrict pD;
    const xtfloat * restrict pA;
    for (n = 0; n<10; n++)
        for (p = 0; p<10; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloat *)(Z + n + p * 10);
            XT_LSIP(Yw, pA, SZ_F32);
            pY = (xtfloat*)(y + p);
            for (m = 0; m<n; m++)
            {
                XT_LSXP(Yr, pY, 10 * SZ_F32);
                XT_LSIP(R0, pR, SZ_F32);
                XT_MSUB_S(Yw, R0, Yr);
            }
            pD = (const xtfloat *)(D + n);
            XT_LSIP(D0, pD, SZ_F32);
            Yw = XT_MUL_S(Yw, D0);
            XT_SSIP(Yw, pY, SZ_F32);
        }
#else
    int n, m;
    xtfloat R0, D0;
    xtfloat Yr0, Yr1, Yr2, Yr3, Yr4, Yr5, Yr6, Yr7, Yr8, Yr9;
    xtfloat Yw0, Yw1, Yw2, Yw3, Yw4, Yw5, Yw6, Yw7, Yw8, Yw9;
    xtfloat * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)R;
    const xtfloat * restrict pD = (const xtfloat *)D;
    const xtfloat * restrict pA = (const xtfloat *)Z;
    for (n = 0; n<10; n++)
    {
        pY = (xtfloat*)y;
        XT_LSIP(Yw0, pA, 10 * SZ_F32);
        XT_LSIP(Yw1, pA, 10 * SZ_F32);
        XT_LSIP(Yw2, pA, 10 * SZ_F32);
        XT_LSIP(Yw3, pA, 10 * SZ_F32);
        XT_LSIP(Yw4, pA, 10 * SZ_F32);
        XT_LSIP(Yw5, pA, 10 * SZ_F32);
        XT_LSIP(Yw6, pA, 10 * SZ_F32);
        XT_LSIP(Yw7, pA, 10 * SZ_F32);
        XT_LSIP(Yw8, pA, 10 * SZ_F32);
        XT_LSXP(Yw9, pA, (-10 * (10 - 1) + 1)*SZ_F32);

        for (m = 0; m<n; m++)
        {
            XT_LSIP(Yr0, pY, SZ_F32);
            XT_LSIP(Yr1, pY, SZ_F32);
            XT_LSIP(Yr2, pY, SZ_F32);
            XT_LSIP(Yr3, pY, SZ_F32);
            XT_LSIP(Yr4, pY, SZ_F32);
            XT_LSIP(Yr5, pY, SZ_F32);
            XT_LSIP(Yr6, pY, SZ_F32);
            XT_LSIP(Yr7, pY, SZ_F32);
            XT_LSIP(Yr8, pY, SZ_F32);
            XT_LSIP(Yr9, pY, SZ_F32);
            XT_LSIP(R0, pR, SZ_F32);
            XT_MSUB_S(Yw0, R0, Yr0);
            XT_MSUB_S(Yw1, R0, Yr1);
            XT_MSUB_S(Yw2, R0, Yr2);
            XT_MSUB_S(Yw3, R0, Yr3);
            XT_MSUB_S(Yw4, R0, Yr4);
            XT_MSUB_S(Yw5, R0, Yr5);
            XT_MSUB_S(Yw6, R0, Yr6);
            XT_MSUB_S(Yw7, R0, Yr7);
            XT_MSUB_S(Yw8, R0, Yr8);
            XT_MSUB_S(Yw9, R0, Yr9);
        }

        XT_LSIP(D0, pD, SZ_F32);
        Yw0 = XT_MUL_S(Yw0, D0);
        Yw1 = XT_MUL_S(Yw1, D0);
        Yw2 = XT_MUL_S(Yw2, D0);
        Yw3 = XT_MUL_S(Yw3, D0);
        Yw4 = XT_MUL_S(Yw4, D0);
        Yw5 = XT_MUL_S(Yw5, D0);
        Yw6 = XT_MUL_S(Yw6, D0);
        Yw7 = XT_MUL_S(Yw7, D0);
        Yw8 = XT_MUL_S(Yw8, D0);
        Yw9 = XT_MUL_S(Yw9, D0);

        XT_SSIP(Yw0, pY, SZ_F32);
        XT_SSIP(Yw1, pY, SZ_F32);
        XT_SSIP(Yw2, pY, SZ_F32);
        XT_SSIP(Yw3, pY, SZ_F32);
        XT_SSIP(Yw4, pY, SZ_F32);
        XT_SSIP(Yw5, pY, SZ_F32);
        XT_SSIP(Yw6, pY, SZ_F32);
        XT_SSIP(Yw7, pY, SZ_F32);
        XT_SSIP(Yw8, pY, SZ_F32);
        XT_SSIP(Yw9, pY, SZ_F32);
        pR += 1;
    }
#endif
}
/*
backward recursion: P!=1
*/
static void rbkw10x10f(
    float32_t* restrict x,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k;
    xtfloat * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 10 * (10 + 1) / 2 - 1); //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (10 - 1));
    const xtfloat * restrict pY = (xtfloat*)(y + (10 * 10 - 1));
    xtfloat Xw0, Xw1, Xw2, Xw3, Xw4, Xw5, Xw6, Xw7, Xw8, Xw9;
    xtfloat Xr0, Xr1, Xr2, Xr3, Xr4, Xr5, Xr6, Xr7, Xr8, Xr9;
    xtfloat R0, D0;

    for (k = 10 - 1; k >= 0; k--)
    {
        {
            pX = (xtfloat*)(x + 10 * 10 - 1); // last element in X
            pR = (xtfloat*)(tpR--); //points to the end of row k in R
            XT_LSXP(Xw9, pY, -SZ_F32);
            XT_LSXP(Xw8, pY, -SZ_F32);
            XT_LSXP(Xw7, pY, -SZ_F32);
            XT_LSXP(Xw6, pY, -SZ_F32);
            XT_LSXP(Xw5, pY, -SZ_F32);
            XT_LSXP(Xw4, pY, -SZ_F32);
            XT_LSXP(Xw3, pY, -SZ_F32);
            XT_LSXP(Xw2, pY, -SZ_F32);
            XT_LSXP(Xw1, pY, -SZ_F32);
            XT_LSXP(Xw0, pY, -SZ_F32);

            for (m = 0; m < 10 - k - 1; m++)
            {
                XT_LSXP(Xr9, pX, -SZ_F32);
                XT_LSXP(Xr8, pX, -SZ_F32);
                XT_LSXP(Xr7, pX, -SZ_F32);
                XT_LSXP(Xr6, pX, -SZ_F32);
                XT_LSXP(Xr5, pX, -SZ_F32);
                XT_LSXP(Xr4, pX, -SZ_F32);
                XT_LSXP(Xr3, pX, -SZ_F32);
                XT_LSXP(Xr2, pX, -SZ_F32);
                XT_LSXP(Xr1, pX, -SZ_F32);
                XT_LSXP(Xr0, pX, -SZ_F32);

                XT_LSXP(R0, pR, -(10 - 1 - m)*SZ_F32);
                XT_MSUB_S(Xw9, Xr9, R0);
                XT_MSUB_S(Xw8, Xr8, R0);
                XT_MSUB_S(Xw7, Xr7, R0);
                XT_MSUB_S(Xw6, Xr6, R0);
                XT_MSUB_S(Xw5, Xr5, R0);
                XT_MSUB_S(Xw4, Xr4, R0);
                XT_MSUB_S(Xw3, Xr3, R0);
                XT_MSUB_S(Xw2, Xr2, R0);
                XT_MSUB_S(Xw1, Xr1, R0);
                XT_MSUB_S(Xw0, Xr0, R0);
            }
            XT_LSXP(D0, pD, -SZ_F32);
            Xw9 = XT_MUL_S(Xw9, D0);
            Xw8 = XT_MUL_S(Xw8, D0);
            Xw7 = XT_MUL_S(Xw7, D0);
            Xw6 = XT_MUL_S(Xw6, D0);
            Xw5 = XT_MUL_S(Xw5, D0);
            Xw4 = XT_MUL_S(Xw4, D0);
            Xw3 = XT_MUL_S(Xw3, D0);
            Xw2 = XT_MUL_S(Xw2, D0);
            Xw1 = XT_MUL_S(Xw1, D0);
            Xw0 = XT_MUL_S(Xw0, D0);

            XT_SSXP(Xw9, pX, -SZ_F32);
            XT_SSXP(Xw8, pX, -SZ_F32);
            XT_SSXP(Xw7, pX, -SZ_F32);
            XT_SSXP(Xw6, pX, -SZ_F32);
            XT_SSXP(Xw5, pX, -SZ_F32);
            XT_SSXP(Xw4, pX, -SZ_F32);
            XT_SSXP(Xw3, pX, -SZ_F32);
            XT_SSXP(Xw2, pX, -SZ_F32);
            XT_SSXP(Xw1, pX, -SZ_F32);
            XT_SSXP(Xw0, pX, -SZ_F32);
        }
    }
}
/*-------------------------------------------------------------------------
Matrix (Pseudo) Inversion
Obtain Left Inverse of a matrix using Cholesky Decomposition
The result is matrix x = A^-1
Fixed point API requires explicit setting of fixed point representation of
input/output matrices as well as for internal temporary matrices such as R
(Cholesky decomposition) and y (decision of R'y=(A'*B))


Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]      matrix A, for fixed point API, the representation is Q(qA)
sigma2      regularization term, for fixed point API, the
representation is Q(2*qA-30)
qRA         qR-qA; difference between fixed point representations of R
and A matrices (for the fixed point API only). Should be
equal or less than 0 (typically -2).
qYBRA       qY-qB+qR-qA, combination of fixed point representations of
matrices y, B, R and A (for the fixed point API only). Since
for matrix inversion we simply use identity matrix B, we may
always suppose qB=31
qXYR        combination of fixed point representation (matrices R, x and y)
qXYR=qX-qY+qR (for the fixed point API only)
Output:
x[N*M]      Left Inverse of the matrix A, for fixed point API, the
representation is Q(qX)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  matcholpseudoinv10x10f(void* pScr,
    float32_t *x,
    const float32_t * A,
    const float32_t  sigma2)
{
    float32_t * restrict D; int SD;
    float32_t * restrict R; int SR;
    float32_t * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);
    R = (float32_t *)pScr;
    D = (float32_t *)(((uintptr_t)R) + SR);
    y = (float32_t *)(((uintptr_t)D) + SD);
    pScr = (float32_t *)(((uintptr_t)y) + SY);
    matcholdecomp10x10f(pScr, R, D, A, sigma2);
    rfwd10x10f(y, R, D, A);
    rbkw10x10f(x, R, D, y);
}

size_t   matcholpseudoinv10x10f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);

    s_dc = matcholdecomp10x10f_getScratchSize();

    return SD + SR + SY + s_dc;
}
#else
DISCARD_FUN(void, matcholpseudoinv10x10f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

  size_t  matcholpseudoinv10x10f_getScratchSize()
{
	return 0;
}
#endif
