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
static void rfwd4x4f(
    float32_t* restrict y,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict Z)
{
#if 0
    int n, m, p;
    xtfloat Yw, Yr, R0, D0, A0;
    xtfloat * restrict pY;
    const xtfloat * restrict pR;
    const xtfloat * restrict pD;
    const xtfloat * restrict pA;
    for (n = 0; n<4; n++)
        for (p = 0; p<4; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1))/2);
            pA = (const xtfloat *)(Z + n + p * 4);
            XT_LSIP(Yw, pA, SZ_F32);
            pY = (xtfloat*)(y + p);
            for (m = 0; m<n; m++)
            {
                XT_LSXP(Yr, pY, 4*SZ_F32);
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
    xtfloatx2 Yr01, Yr23;
    xtfloat Yw0, Yw1, Yw2, Yw3;
    xtfloatx2 Yw01, Yw23;
    xtfloatx2 * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)R;
    const xtfloat * restrict pD = (const xtfloat *)D;
    const xtfloat * restrict pA = (const xtfloat *)Z;
    for (n = 0; n<4; n++)
    {
        pY = (xtfloatx2*)y;
        XT_LSIP(Yw0, pA, 4 * SZ_F32);
        XT_LSIP(Yw1, pA, 4 * SZ_F32);
        XT_LSIP(Yw2, pA, 4 * SZ_F32);
        XT_LSXP(Yw3, pA, (-4 * 3 + 1)*SZ_F32);
        Yw01 = XT_SEL32_HH_SX2(Yw0, Yw1);
        Yw23 = XT_SEL32_HH_SX2(Yw2, Yw3);

        for (m = 0; m<n; m++)
        {
            XT_LSX2IP(Yr01, pY, 2 * SZ_F32);
            XT_LSX2IP(Yr23, pY, 2 * SZ_F32);
            XT_LSIP(R0, pR, SZ_F32);

            XT_MSUB_SX2(Yw01, R0, Yr01);
            XT_MSUB_SX2(Yw23, R0, Yr23);
        }

        XT_LSIP(D0, pD, SZ_F32);
        Yw01 = XT_MUL_SX2(Yw01, D0);
        Yw23 = XT_MUL_SX2(Yw23, D0);

        XT_SSX2IP(Yw01, pY, 2 * SZ_F32);
        XT_SSX2IP(Yw23, pY, 2 * SZ_F32);
        pR += 1;
    }
#endif
}
/*
backward recursion: P!=1
*/
#if 0 // basic, reqire transformed R
static void rbkw4x4f(
    float32_t* restrict x,
    const float32_t* restrict Rt,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k, p;
    const xtfloat* restrict pRt;
    xtfloat * restrict pXw;
    xtfloat * restrict pXr;
    xtfloat * restrict pD = (xtfloat*)(D+(4-1));
    xtfloat * restrict pY;

    xtfloat X, Xw, R, D0;

    for (k = 4 - 1; k >= 0; k--)
    {
        pY = (xtfloat*)(y + k*4);
        pXw = (xtfloat*)(x + k*4);
        for (p = 0; p < 4; p++)
        {

            XT_LSIP(Xw, pY, SZ_F32);

            pRt = (xtfloat*)(Rt + (4 - k - 1)*(4 - k - 2)/2);
            for (m = 0; m < 4 - k - 1; m++)
            {
                pXr = (xtfloat*)(x + (k + 1 + m)*4 + p);
                XT_LSIP(X, pXr, 0);
                XT_LSIP(R, pRt, SZ_F32);
                XT_MSUB_S(Xw, X, R);
            }
            XT_LSIP(D0, pD, 0);
            Xw = XT_MUL_S(Xw, D0);

            XT_SSIP(Xw, pXw, SZ_F32);
        }
        pD--;
    }
}
#else // does not require tranformed R, reverse inner loop
static void rbkw4x4f(
    float32_t* restrict x,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k;
    xtfloatx2 * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 4 * 5 / 2 - 1); //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (4 - 1));
    const xtfloatx2 * restrict pY = (xtfloatx2*)(y + (4 * 4 - 2));
    xtfloatx2 Xw01, Xw23;
    xtfloatx2 Xr01, Xr23;
    xtfloat R0, D0;

    for (k = 4 - 1; k >= 0; k--)
    {
        {
            pX = (xtfloatx2*)(x + 4 * 4 - 2); // last 2 elements in X
            pR = (xtfloat*)(tpR--); //points to the end of row k in R
            XT_LSX2XP(Xw23, pY, -2 * SZ_F32);
            XT_LSX2XP(Xw01, pY, -2 * SZ_F32);

            for (m = 0; m < 4 - k - 1; m++)
            {
                XT_LSX2XP(Xr23, pX, -2 * SZ_F32);
                XT_LSX2XP(Xr01, pX, -2 * SZ_F32);

                XT_LSXP(R0, pR, -(3 - m)*SZ_F32);
                XT_MSUB_SX2(Xw23, Xr23, R0);
                XT_MSUB_SX2(Xw01, Xr01, R0);
            }
            XT_LSXP(D0, pD, -SZ_F32);
            Xw23 = XT_MUL_SX2(Xw23, D0);
            Xw01 = XT_MUL_SX2(Xw01, D0);
            XT_SSX2XP(Xw23, pX, -2 * SZ_F32);
            XT_SSX2XP(Xw01, pX, -2 * SZ_F32);

        }
    }
#endif
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
void  matcholpseudoinv4x4f(void* pScr,
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

    SD = 4;
    SR = (((4 + 1)*4) >> 1);
    SY = 4*4;
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
    matcholdecomp4x4f(pScr, R, D, A, sigma2);
    rfwd4x4f(y, R, D, A);
    rbkw4x4f(x, R, D, y);
}

size_t   matcholpseudoinv4x4f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 4;
    SR = (((4 + 1)*4) >> 1);
    SY = 4*4;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);

    s_dc = matcholdecomp4x4f_getScratchSize();

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
static void rfwd4x4f(
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
    for (n = 0; n<4; n++)
        for (p = 0; p<4; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloat *)(Z + n + p * 4);
            XT_LSIP(Yw, pA, SZ_F32);
            pY = (xtfloat*)(y + p);
            for (m = 0; m<n; m++)
            {
                XT_LSXP(Yr, pY, 4 * SZ_F32);
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
    xtfloat Yr0, Yr1, Yr2, Yr3;
    xtfloat Yw0, Yw1, Yw2, Yw3;
    xtfloat * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)R;
    const xtfloat * restrict pD = (const xtfloat *)D;
    const xtfloat * restrict pA = (const xtfloat *)Z;
    for (n = 0; n<4; n++)
    {
        pY = (xtfloat*)y;
        XT_LSIP(Yw0, pA, 4 * SZ_F32);
        XT_LSIP(Yw1, pA, 4 * SZ_F32);
        XT_LSIP(Yw2, pA, 4 * SZ_F32);
        XT_LSXP(Yw3, pA, (-4 * 3 + 1)*SZ_F32);

        for (m = 0; m<n; m++)
        {
            XT_LSIP(Yr0, pY, SZ_F32);
            XT_LSIP(Yr1, pY, SZ_F32);
            XT_LSIP(Yr2, pY, SZ_F32);
            XT_LSIP(Yr3, pY, SZ_F32);
            XT_LSIP(R0, pR, SZ_F32);
            XT_MSUB_S(Yw0, R0, Yr0);
            XT_MSUB_S(Yw1, R0, Yr1);
            XT_MSUB_S(Yw2, R0, Yr2);
            XT_MSUB_S(Yw3, R0, Yr3);
        }

        XT_LSIP(D0, pD, SZ_F32);
        Yw0 = XT_MUL_S(Yw0, D0);
        Yw1 = XT_MUL_S(Yw1, D0);
        Yw2 = XT_MUL_S(Yw2, D0);
        Yw3 = XT_MUL_S(Yw3, D0);

        XT_SSIP(Yw0, pY, SZ_F32);
        XT_SSIP(Yw1, pY, SZ_F32);
        XT_SSIP(Yw2, pY, SZ_F32);
        XT_SSIP(Yw3, pY, SZ_F32);
        pR += 1;
    }
#endif
}
/*
backward recursion: P!=1
*/
#if 0 // basic, reqire transformed R
static void rbkw4x4f(
    float32_t* restrict x,
    const float32_t* restrict Rt,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k, p;
    const xtfloat* restrict pRt;
    xtfloat * restrict pXw;
    xtfloat * restrict pXr;
    xtfloat * restrict pD = (xtfloat*)(D + (4 - 1));
    xtfloat * restrict pY;

    xtfloat X, Xw, R, D0;

    for (k = 4 - 1; k >= 0; k--)
    {
        pY = (xtfloat*)(y + k * 4);
        pXw = (xtfloat*)(x + k * 4);
        for (p = 0; p < 4; p++)
        {

            XT_LSIP(Xw, pY, SZ_F32);

            pRt = (xtfloat*)(Rt + (4 - k - 1)*(4 - k - 2) / 2);
            for (m = 0; m < 4 - k - 1; m++)
            {
                pXr = (xtfloat*)(x + (k + 1 + m) * 4 + p);
                XT_LSIP(X, pXr, 0);
                XT_LSIP(R, pRt, SZ_F32);
                XT_MSUB_S(Xw, X, R);
            }
            XT_LSIP(D0, pD, 0);
            Xw = XT_MUL_S(Xw, D0);

            XT_SSIP(Xw, pXw, SZ_F32);
        }
        pD--;
    }
}
#else // does not require tranformed R, reverse inner loop
static void rbkw4x4f(
    float32_t* restrict x,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k;
    xtfloat * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 4 * 5 / 2 - 1); //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (4 - 1));
    const xtfloat * restrict pY = (xtfloat*)(y + (4 * 4 - 1));
    xtfloat Xw0, Xw1, Xw2, Xw3;
    xtfloat Xr0, Xr1, Xr2, Xr3;
    xtfloat R0, D0;

    for (k = 4 - 1; k >= 0; k--)
    {
        {
            pX = (xtfloat*)(x + 4 * 4 - 1); // last element in X
            pR = (xtfloat*)(tpR--); //points to the end of row k in R
            XT_LSXP(Xw3, pY, -SZ_F32);
            XT_LSXP(Xw2, pY, -SZ_F32);
            XT_LSXP(Xw1, pY, -SZ_F32);
            XT_LSXP(Xw0, pY, -SZ_F32);

            for (m = 0; m < 4 - k - 1; m++)
            {
                XT_LSXP(Xr3, pX, -SZ_F32);
                XT_LSXP(Xr2, pX, -SZ_F32);
                XT_LSXP(Xr1, pX, -SZ_F32);
                XT_LSXP(Xr0, pX, -SZ_F32);

                XT_LSXP(R0, pR, -(3 - m)*SZ_F32);
                XT_MSUB_S(Xw3, Xr3, R0);
                XT_MSUB_S(Xw2, Xr2, R0);
                XT_MSUB_S(Xw1, Xr1, R0);
                XT_MSUB_S(Xw0, Xr0, R0);
            }
            XT_LSXP(D0, pD, -SZ_F32);
            Xw3 = XT_MUL_S(Xw3, D0);
            Xw2 = XT_MUL_S(Xw2, D0);
            Xw1 = XT_MUL_S(Xw1, D0);
            Xw0 = XT_MUL_S(Xw0, D0);

            XT_SSXP(Xw3, pX, -SZ_F32);
            XT_SSXP(Xw2, pX, -SZ_F32);
            XT_SSXP(Xw1, pX, -SZ_F32);
            XT_SSXP(Xw0, pX, -SZ_F32);
        }
    }
#endif
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
void  matcholpseudoinv4x4f(void* pScr,
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

    SD = 4;
    SR = (((4 + 1) * 4) >> 1);
    SY = 4 * 4;
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
    matcholdecomp4x4f(pScr, R, D, A, sigma2);
    rfwd4x4f(y, R, D, A);
    rbkw4x4f(x, R, D, y);
}

size_t   matcholpseudoinv4x4f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 4;
    SR = (((4 + 1) * 4) >> 1);
    SY = 4 * 4;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);

    s_dc = matcholdecomp4x4f_getScratchSize();

    return SD + SR + SY + s_dc;
}

#else
DISCARD_FUN(void, matcholpseudoinv4x4f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

  size_t  matcholpseudoinv4x4f_getScratchSize()
{
	return 0;
}
#endif
