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
#define CF32_WIDTH (HIFI_SIMD_WIDTH/SZ_CF32)
/*
Forward recursion for pseudo inversion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
*/
static void fwd4x4f(
    complex_float* restrict y,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict Z)
{
#if 0
    int n, m, p;
    xtfloatx2 Yw, Yr, R0, D0;
    xtfloatx2 * restrict pY;
    const xtfloatx2 * restrict pR;
    const xtfloatx2 * restrict pD;
    const xtfloatx2 * restrict pA;
    for (n = 0; n<4; n++)
        for (p = 0; p<4; p++)
        {
            pR = (const xtfloatx2 *)(R + (n*(n + 1))/2);
            pA = (const xtfloatx2 *)(Z + n + p * 4);
            XT_LSX2IP(Yw, pA, SZ_CF32);
            Yw = XT_CONJC_S(Yw);
            pY = (xtfloatx2*)(y + p);
            for (m = 0; m<n; m++)
            {
                XT_LSX2XP(Yr, pY, 4*SZ_CF32);
                XT_LSX2IP(R0, pR, SZ_CF32);
                XT_MADDMUX_S(Yw, R0, Yr, 2);
                XT_MADDMUX_S(Yw, R0, Yr, 1);
            }
            pD = (const xtfloatx2 *)(D + n);
            XT_LSX2IP(D0, pD, SZ_CF32);
            Yw = XT_MUL_SX2(Yw, D0);
            XT_SSX2IP(Yw, pY, SZ_CF32);
        }
#else
    int n, m;
    xtfloatx2 R0, D0;
    xtfloatx2 Yr0, Yr1, Yr2, Yr3;
    xtfloatx2 Yw0, Yw1, Yw2, Yw3;
    xtfloatx2 Yw0t, Yw1t, Yw2t, Yw3t;
    xtfloatx2 * restrict pY;
    const xtfloatx2 * restrict pR = (const xtfloatx2 *)R;
    const xtfloatx2 * restrict pD = (const xtfloatx2 *)D;
    const xtfloatx2 * restrict pA = (const xtfloatx2 *)Z;
    for (n = 0; n<4; n++)
        //for (p = 0; p<P; p++)
        {         
            pY = (xtfloatx2*)y;
            XT_LSX2IP(Yw0, pA, 4*SZ_CF32);
            XT_LSX2IP(Yw1, pA, 4*SZ_CF32);
            XT_LSX2IP(Yw2, pA, 4*SZ_CF32);
            XT_LSX2XP(Yw3, pA, (-4*3+1)*SZ_CF32);
            Yw0 = XT_CONJC_S(Yw0);
            Yw1 = XT_CONJC_S(Yw1);
            Yw2 = XT_CONJC_S(Yw2);
            Yw3 = XT_CONJC_S(Yw3);

            Yw0t = Yw1t = Yw2t = Yw3t = XT_CONST_S(0);
            for (m = 0; m<n; m++)
            {
                XT_LSX2IP(Yr0, pY, SZ_CF32);
                XT_LSX2IP(Yr1, pY, SZ_CF32);
                XT_LSX2IP(Yr2, pY, SZ_CF32);
                XT_LSX2IP(Yr3, pY, SZ_CF32);
                XT_LSX2IP(R0, pR, SZ_CF32);

                XT_MADDMUX_S(Yw0, R0, Yr0, 2);
                XT_MADDMUX_S(Yw0t, R0, Yr0, 1);
                XT_MADDMUX_S(Yw1, R0, Yr1, 2);
                XT_MADDMUX_S(Yw1t, R0, Yr1, 1);
                XT_MADDMUX_S(Yw2, R0, Yr2, 2);
                XT_MADDMUX_S(Yw2t, R0, Yr2, 1);
                XT_MADDMUX_S(Yw3, R0, Yr3, 2);
                XT_MADDMUX_S(Yw3t, R0, Yr3, 1);
            }
            Yw0 = XT_ADD_SX2(Yw0, Yw0t);
            Yw1 = XT_ADD_SX2(Yw1, Yw1t);
            Yw2 = XT_ADD_SX2(Yw2, Yw2t);
            Yw3 = XT_ADD_SX2(Yw3, Yw3t);

            XT_LSX2IP(D0, pD, SZ_CF32);
            Yw0 = XT_MUL_SX2(Yw0, D0);
            Yw1 = XT_MUL_SX2(Yw1, D0);
            Yw2 = XT_MUL_SX2(Yw2, D0);
            Yw3 = XT_MUL_SX2(Yw3, D0);
            XT_SSX2IP(Yw0, pY, SZ_CF32);
            XT_SSX2IP(Yw1, pY, SZ_CF32);
            XT_SSX2IP(Yw2, pY, SZ_CF32);
            XT_SSX2IP(Yw3, pY, SZ_CF32);
            pR += 1;
    }
#endif
}
#if 0
/* reversing R matrices for easier readings by rows (diagonal elements are omitted):
original R    transformed R
0 1 3 6 a     d 8 c 4 7 b 1 3 6 a
2 4 7 b
5 8 c
9 d
e

Input:
R[L][SR]        L input matrices
Rt[L*N*(N-1)]   stream of L trasposed matrices
*/
static void transformRf(float32_t* Rt, const float32_t* R, int N)
{
    int n, m;
    const float32_t* pR;
    for (n = 0; n < N; n++)
    {
        pR = R + (N - n)*(N - n + 3) - 2;
        for (m = 0; m < n; m++)
        {
            Rt[0] = pR[0];
            Rt[1] = pR[1];
            pR += 2 * (N - n + m + 1);
            Rt += 2;
        }
    }
}
#endif
/*
backward recursion: P!=1
*/
#if 0
static void bkw4x4f(
    complex_float* restrict x,
    const complex_float* restrict Rt,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k, p;
    const xtfloatx2* restrict pRt;
    xtfloatx2 * restrict pXw;
    xtfloatx2 * restrict pXr;
    xtfloatx2 * restrict pD = (xtfloatx2*)(D+(4-1));
    xtfloatx2 * restrict pY;

    xtfloatx2 X, Xw, R, D0;

    for (k = 4 - 1; k >= 0; k--)
    {
        pY = (xtfloatx2*)(y + k*4);
        pXw = (xtfloatx2*)(x + k*4);
        // calculate y(m,:)-R(m,:)*X, 1xP
        for (p = 0; p < 4; p++)
        {
            
            XT_LSX2IP(Xw, pY, SZ_CF32);
            
            pRt = (xtfloatx2*)(Rt + (4 - k - 1)*(4 - k - 2)/2);
            for (m = 0; m < 4 - k - 1; m++)
            {
                pXr = (xtfloatx2*)(x + (k + 1 + m)*4 + p);
                XT_LSX2IP(X, pXr, 0);
                XT_LSX2IP(R, pRt, SZ_CF32);
                XT_MADDMUX_S(Xw, X, R, 2);
                XT_MADDMUX_S(Xw, X, R, 3);
            }
            //pD = (xtfloatx2*)(D + k);
            XT_LSX2IP(D0, pD, 0);
            Xw = XT_MUL_SX2(Xw, D0);
            
            XT_SSX2IP(Xw, pXw, SZ_CF32);
        }
        pD--;
    }
}
#elif 0 // require transformed R
static void bkw4x4f(
    complex_float* restrict x,
    const complex_float* restrict Rt,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k;
    const xtfloatx2* restrict pRt = (xtfloatx2*)Rt;
    xtfloatx2 * restrict pXw = (xtfloatx2*)(x + (4*4 - 1));
    xtfloatx2 * restrict pXr;
    xtfloatx2 * restrict pD = (xtfloatx2*)(D+(4-1));
    xtfloatx2 * restrict pY = (xtfloatx2*)(y + (4*4 - 1));
    xtfloatx2 X0, X1, X2, X3;
    xtfloatx2 Xw0, Xw1, Xw2, Xw3;
    xtfloatx2 Xw0t, Xw1t, Xw2t, Xw3t;
    xtfloatx2 R0, D0;

    for (k = 4 - 1; k >= 0; k--)
    {
        //pY = (xtfloatx2*)(y + 2 * k*P + 2*(P-1));
        //pXw = (xtfloatx2*)(x + 2 * k*P + 2*(P-1));
        //for (p = 0; p < P; p++)
        {
            XT_LSX2XP(Xw3, pY, -SZ_CF32);
            XT_LSX2XP(Xw2, pY, -SZ_CF32);
            XT_LSX2XP(Xw1, pY, -SZ_CF32);
            XT_LSX2XP(Xw0, pY, -SZ_CF32);

            pXr = (xtfloatx2*)(x + (k + 1)*4);
            //pXr = (xtfloatx2*)(x + (k + 1)*4 + 4*(4 - k - 1) - 1);
            //pRt = (xtfloatx2*)(Rt + (4 - k - 1)*(4 - k - 2) / 2 + (4 - k - 1) - 1);
            Xw0t = Xw1t = Xw2t = Xw3t = XT_CONST_S(0);
            for (m = 0; m < 4 - k - 1; m++)
            {
                XT_LSX2XP(X0, pXr, SZ_CF32);
                XT_LSX2XP(X1, pXr, SZ_CF32);
                XT_LSX2XP(X2, pXr, SZ_CF32);
                XT_LSX2XP(X3, pXr, SZ_CF32);
                XT_LSX2XP(R0, pRt, SZ_CF32);
                XT_MADDMUX_S(Xw0, X0, R0, 2);
                XT_MADDMUX_S(Xw0t, X0, R0, 3);
                XT_MADDMUX_S(Xw1, X1, R0, 2);
                XT_MADDMUX_S(Xw1t, X1, R0, 3);
                XT_MADDMUX_S(Xw2, X2, R0, 2);
                XT_MADDMUX_S(Xw2t, X2, R0, 3);
                XT_MADDMUX_S(Xw3, X3, R0, 2);
                XT_MADDMUX_S(Xw3t, X3, R0, 3);
            }
            Xw0 = XT_ADD_SX2(Xw0, Xw0t);
            Xw1 = XT_ADD_SX2(Xw1, Xw1t);
            Xw2 = XT_ADD_SX2(Xw2, Xw2t);
            Xw3 = XT_ADD_SX2(Xw3, Xw3t);

            XT_LSX2XP(D0, pD, -SZ_CF32);
            Xw0 = XT_MUL_SX2(Xw0, D0);
            Xw1 = XT_MUL_SX2(Xw1, D0);
            Xw2 = XT_MUL_SX2(Xw2, D0);
            Xw3 = XT_MUL_SX2(Xw3, D0);

            XT_SSX2XP(Xw3, pXw, -SZ_CF32);
            XT_SSX2XP(Xw2, pXw, -SZ_CF32);
            XT_SSX2XP(Xw1, pXw, -SZ_CF32);
            XT_SSX2XP(Xw0, pXw, -SZ_CF32);
        }
    }
}
#else // does not require tranformed R, reverse inner loop
static void bkw4x4f(
    complex_float* restrict x,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k;
    xtfloatx2 * restrict pX;
    const xtfloatx2 * restrict pR;
    const xtfloatx2 * restrict tpR = (xtfloatx2*)(R + 4 * 5 / 2 - 1); //last element in R
    const xtfloatx2 * restrict pD = (xtfloatx2*)(D + (4 - 1));
    const xtfloatx2 * restrict pY = (xtfloatx2*)(y + (4 * 4 - 1));
    xtfloatx2 X0, X1, X2, X3;
    xtfloatx2 Xw0, Xw1, Xw2, Xw3;
    xtfloatx2 Xw0t, Xw1t, Xw2t, Xw3t;
    xtfloatx2 R0, D0;
    
    for (k = 4 - 1; k >= 0; k--)
    {
        //for (p = 0; p < P; p++)
        {
            pX = (xtfloatx2*)(x + 4 * 4 - 1); // last element in X
            pR = (xtfloatx2*)(tpR--); //points to the end of row k in R
            XT_LSX2XP(Xw3, pY, -SZ_CF32);
            XT_LSX2XP(Xw2, pY, -SZ_CF32);
            XT_LSX2XP(Xw1, pY, -SZ_CF32);
            XT_LSX2XP(Xw0, pY, -SZ_CF32);
            
            Xw0t = Xw1t = Xw2t = Xw3t = XT_CONST_S(0);
            for (m = 0; m < 4 - k - 1; m++)
            {
                XT_LSX2XP(X3, pX, -SZ_CF32);
                XT_LSX2XP(X2, pX, -SZ_CF32);
                XT_LSX2XP(X1, pX, -SZ_CF32);
                XT_LSX2XP(X0, pX, -SZ_CF32);
                XT_LSX2XP(R0, pR, -(3 - m)*SZ_CF32);
                XT_MADDMUX_S(Xw0, X0, R0, 2);
                XT_MADDMUX_S(Xw0t, X0, R0, 3);
                XT_MADDMUX_S(Xw1, X1, R0, 2);
                XT_MADDMUX_S(Xw1t, X1, R0, 3);
                XT_MADDMUX_S(Xw2, X2, R0, 2);
                XT_MADDMUX_S(Xw2t, X2, R0, 3);
                XT_MADDMUX_S(Xw3, X3, R0, 2);
                XT_MADDMUX_S(Xw3t, X3, R0, 3);
            }
            Xw0 = XT_ADD_SX2(Xw0, Xw0t);
            Xw1 = XT_ADD_SX2(Xw1, Xw1t);
            Xw2 = XT_ADD_SX2(Xw2, Xw2t);
            Xw3 = XT_ADD_SX2(Xw3, Xw3t);

            XT_LSX2XP(D0, pD, -SZ_CF32);
            Xw0 = XT_MUL_SX2(Xw0, D0);
            Xw1 = XT_MUL_SX2(Xw1, D0);
            Xw2 = XT_MUL_SX2(Xw2, D0);
            Xw3 = XT_MUL_SX2(Xw3, D0);

            XT_SSX2XP(Xw3, pX, -SZ_CF32);
            XT_SSX2XP(Xw2, pX, -SZ_CF32);
            XT_SSX2XP(Xw1, pX, -SZ_CF32);
            XT_SSX2XP(Xw0, pX, -SZ_CF32);
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
void  cmatcholpseudoinv4x4f(void* pScr,
    complex_float *x,
    const complex_float * A,
    const float32_t  sigma2)
{
    complex_float * restrict D; int SD;
    complex_float * restrict R; int SR;
    complex_float * restrict y; int SY;

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
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);
    R = (complex_float *)pScr;
    D = (complex_float *)(((uintptr_t)R) + SR);
    y = (complex_float *)(((uintptr_t)D) + SD);
    pScr = (complex_float *)(((uintptr_t)y) + SY);
    
    cmatcholdecomp4x4f(pScr, R, D, A, sigma2);
    fwd4x4f(y, R, D, A);
    //transformRf((float32_t*)Rt, (float32_t*)R, 4);
    bkw4x4f(x, R, D, y);
}

size_t  cmatcholpseudoinv4x4f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY/*, SRt*/;

    SD = 4;
    SR = (((4 + 1)*4) >> 1);
    SY = 4*4;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);

    s_dc = cmatcholdecomp4x4f_getScratchSize();
    return SD + SR + SY + s_dc;
}
#elif (HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
#define CF32_WIDTH (HIFI_SIMD_WIDTH/SZ_CF32)
/*
Forward recursion for pseudo inversion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
*/
static void fwd4x4f(
    complex_float* restrict y,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict Z)
{
#if 0
    int n, m, p;
    xtfloat Yw_re, Yw_im, Yr_re, Yr_im, R_re, R_im, D_re, D_im;
    xtfloat * restrict pY;
    const xtfloat * restrict pR;
    const xtfloat * restrict pD;
    const xtfloat * restrict pA;
    for (n = 0; n<4; n++)
    {
        //pA = (const xtfloat *)(Z + n);
        for (p = 0; p < 4; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloat *)(Z + n + p * 4);
            XT_LSIP(Yw_re, pA, SZ_F32);
            XT_LSIP(Yw_im, pA, 4 * SZ_CF32 - SZ_F32);
            Yw_im = XT_NEG_S(Yw_im);
            pY = (xtfloat*)(y + p);
            for (m = 0; m < n; m++)
            {
                XT_LSXP(Yr_re, pY, SZ_F32);
                XT_LSXP(Yr_im, pY, 4 * SZ_CF32 - SZ_F32);
                XT_LSIP(R_re, pR, SZ_F32);
                XT_LSIP(R_im, pR, SZ_F32);
                XT_MSUB_S(Yw_re, R_re, Yr_re);
                XT_MSUB_S(Yw_re, R_im, Yr_im);
                XT_MSUB_S(Yw_im, R_re, Yr_im);
                XT_MADD_S(Yw_im, R_im, Yr_re);

                //XT_MADDMUX_S(Yw, R0, Yr, 2);
                //XT_MADDMUX_S(Yw, R0, Yr, 1);
            }
            pD = (const xtfloat *)(D + n);
            XT_LSIP(D_re, pD, SZ_F32);
            XT_LSIP(D_im, pD, SZ_F32);
            Yw_re = XT_MUL_S(Yw_re, D_re);
            Yw_im = XT_MUL_S(Yw_im, D_im);
            XT_SSIP(Yw_re, pY, SZ_F32);
            XT_SSIP(Yw_im, pY, SZ_F32);
        }
    }
#elif 1
    int n, m;
    xtfloat R_re, R_im, D_re, D_im;
    xtfloat Yw0_re, Yw1_re, Yw2_re, Yw3_re;
    xtfloat Yw0_im, Yw1_im, Yw2_im, Yw3_im;
    xtfloat Yr0_re, Yr1_re, Yr2_re, Yr3_re;
    xtfloat Yr0_im, Yr1_im, Yr2_im, Yr3_im;

    xtfloat * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)(R);
    const xtfloat * restrict pD = (const xtfloat *)(D);
    const xtfloat * restrict pA;
    for (n = 0; n<4; n++)
    {
        pA = (const xtfloat *)(Z + n);
        //pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
        XT_LSIP(Yw0_re, pA, SZ_F32);
        XT_LSIP(Yw0_im, pA, 4 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw1_re, pA, SZ_F32);
        XT_LSIP(Yw1_im, pA, 4 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw2_re, pA, SZ_F32);
        XT_LSIP(Yw2_im, pA, 4 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw3_re, pA, SZ_F32);
        XT_LSIP(Yw3_im, pA, 4 * SZ_CF32 - SZ_F32);
        Yw0_im = XT_NEG_S(Yw0_im);
        Yw1_im = XT_NEG_S(Yw1_im);
        Yw2_im = XT_NEG_S(Yw2_im);
        Yw3_im = XT_NEG_S(Yw3_im);
        pY = (xtfloat*)(y);
        for (m = 0; m<n; m++)
        {
            XT_LSXP(Yr0_re, pY, SZ_F32);
            XT_LSXP(Yr0_im, pY, SZ_F32);
            XT_LSXP(Yr1_re, pY, SZ_F32);
            XT_LSXP(Yr1_im, pY, SZ_F32);
            XT_LSXP(Yr2_re, pY, SZ_F32);
            XT_LSXP(Yr2_im, pY, SZ_F32);
            XT_LSXP(Yr3_re, pY, SZ_F32);
            XT_LSXP(Yr3_im, pY, SZ_F32);
            XT_LSIP(R_re, pR, SZ_F32);
            XT_LSIP(R_im, pR, SZ_F32);
            XT_MSUB_S(Yw0_re, R_re, Yr0_re);
            XT_MSUB_S(Yw0_re, R_im, Yr0_im);
            XT_MSUB_S(Yw0_im, R_re, Yr0_im);
            XT_MADD_S(Yw0_im, R_im, Yr0_re);
            XT_MSUB_S(Yw1_re, R_re, Yr1_re);
            XT_MSUB_S(Yw1_re, R_im, Yr1_im);
            XT_MSUB_S(Yw1_im, R_re, Yr1_im);
            XT_MADD_S(Yw1_im, R_im, Yr1_re);
            XT_MSUB_S(Yw2_re, R_re, Yr2_re);
            XT_MSUB_S(Yw2_re, R_im, Yr2_im);
            XT_MSUB_S(Yw2_im, R_re, Yr2_im);
            XT_MADD_S(Yw2_im, R_im, Yr2_re);
            XT_MSUB_S(Yw3_re, R_re, Yr3_re);
            XT_MSUB_S(Yw3_re, R_im, Yr3_im);
            XT_MSUB_S(Yw3_im, R_re, Yr3_im);
            XT_MADD_S(Yw3_im, R_im, Yr3_re);
        }
        //pD = (const xtfloat *)(D + n);
        XT_LSIP(D_re, pD, SZ_F32);
        XT_LSIP(D_im, pD, SZ_F32);
        Yw0_re = XT_MUL_S(Yw0_re, D_re);
        Yw0_im = XT_MUL_S(Yw0_im, D_im);
        Yw1_re = XT_MUL_S(Yw1_re, D_re);
        Yw1_im = XT_MUL_S(Yw1_im, D_im);
        Yw2_re = XT_MUL_S(Yw2_re, D_re);
        Yw2_im = XT_MUL_S(Yw2_im, D_im);
        Yw3_re = XT_MUL_S(Yw3_re, D_re);
        Yw3_im = XT_MUL_S(Yw3_im, D_im);
        XT_SSIP(Yw0_re, pY, SZ_F32);
        XT_SSIP(Yw0_im, pY, SZ_F32);
        XT_SSIP(Yw1_re, pY, SZ_F32);
        XT_SSIP(Yw1_im, pY, SZ_F32);
        XT_SSIP(Yw2_re, pY, SZ_F32);
        XT_SSIP(Yw2_im, pY, SZ_F32);
        XT_SSIP(Yw3_re, pY, SZ_F32);
        XT_SSIP(Yw3_im, pY, SZ_F32);
        pR += 2;
    }
#endif
}
/*
backward recursion: P!=1
does not require tranformed R, reverse inner loop
*/
static void bkw4x4f(
    complex_float* restrict x,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k;
    xtfloat * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 4 * 5 / 2 - 1) + 1; //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (4 - 1)) + 1;
    const xtfloat * restrict pY = (xtfloat*)(y + (4 * 4 - 1)) + 1;
    xtfloat X0_re, X1_re, X2_re, X3_re;
    xtfloat X0_im, X1_im, X2_im, X3_im;
    xtfloat Xw0_re, Xw1_re, Xw2_re, Xw3_re;
    xtfloat Xw0_im, Xw1_im, Xw2_im, Xw3_im;
    xtfloat R0_re, R0_im, D0_re, D0_im;

    for (k = 4 - 1; k >= 0; k--)
    {
        //for (p = 0; p < P; p++)
        {
            pX = (xtfloat*)(x + 4 * 4 - 1) + 1; // last element in X
            pR = (xtfloat*)tpR; //points to the end of row k in R
            tpR -= 2;
            XT_LSXP(Xw3_im, pY, -SZ_F32);
            XT_LSXP(Xw3_re, pY, -SZ_F32);
            XT_LSXP(Xw2_im, pY, -SZ_F32);
            XT_LSXP(Xw2_re, pY, -SZ_F32);
            XT_LSXP(Xw1_im, pY, -SZ_F32);
            XT_LSXP(Xw1_re, pY, -SZ_F32);
            XT_LSXP(Xw0_im, pY, -SZ_F32);
            XT_LSXP(Xw0_re, pY, -SZ_F32);

            for (m = 0; m < 4 - k - 1; m++)
            {
                XT_LSXP(X3_im, pX, -SZ_F32);
                XT_LSXP(X3_re, pX, -SZ_F32);
                XT_LSXP(X2_im, pX, -SZ_F32);
                XT_LSXP(X2_re, pX, -SZ_F32);
                XT_LSXP(X1_im, pX, -SZ_F32);
                XT_LSXP(X1_re, pX, -SZ_F32);
                XT_LSXP(X0_im, pX, -SZ_F32);
                XT_LSXP(X0_re, pX, -SZ_F32);
                XT_LSXP(R0_im, pR, -SZ_F32);
                XT_LSXP(R0_re, pR, -(3 - m)*SZ_CF32 + SZ_F32);
                XT_MSUB_S(Xw0_re, X0_re, R0_re);
                XT_MADD_S(Xw0_re, X0_im, R0_im);
                XT_MSUB_S(Xw0_im, X0_re, R0_im);
                XT_MSUB_S(Xw0_im, X0_im, R0_re);
                XT_MSUB_S(Xw1_re, X1_re, R0_re);
                XT_MADD_S(Xw1_re, X1_im, R0_im);
                XT_MSUB_S(Xw1_im, X1_re, R0_im);
                XT_MSUB_S(Xw1_im, X1_im, R0_re);
                XT_MSUB_S(Xw2_re, X2_re, R0_re);
                XT_MADD_S(Xw2_re, X2_im, R0_im);
                XT_MSUB_S(Xw2_im, X2_re, R0_im);
                XT_MSUB_S(Xw2_im, X2_im, R0_re);
                XT_MSUB_S(Xw3_re, X3_re, R0_re);
                XT_MADD_S(Xw3_re, X3_im, R0_im);
                XT_MSUB_S(Xw3_im, X3_re, R0_im);
                XT_MSUB_S(Xw3_im, X3_im, R0_re);
            }
            XT_LSXP(D0_im, pD, -SZ_F32);
            XT_LSXP(D0_re, pD, -SZ_F32);
            Xw0_re = XT_MUL_S(Xw0_re, D0_re);
            Xw0_im = XT_MUL_S(Xw0_im, D0_im);
            Xw1_re = XT_MUL_S(Xw1_re, D0_re);
            Xw1_im = XT_MUL_S(Xw1_im, D0_im);
            Xw2_re = XT_MUL_S(Xw2_re, D0_re);
            Xw2_im = XT_MUL_S(Xw2_im, D0_im);
            Xw3_re = XT_MUL_S(Xw3_re, D0_re);
            Xw3_im = XT_MUL_S(Xw3_im, D0_im);
            
            XT_SSXP(Xw3_im, pX, -SZ_F32);
            XT_SSXP(Xw3_re, pX, -SZ_F32);
            XT_SSXP(Xw2_im, pX, -SZ_F32);
            XT_SSXP(Xw2_re, pX, -SZ_F32);
            XT_SSXP(Xw1_im, pX, -SZ_F32);
            XT_SSXP(Xw1_re, pX, -SZ_F32);
            XT_SSXP(Xw0_im, pX, -SZ_F32);
            XT_SSXP(Xw0_re, pX, -SZ_F32);
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
void  cmatcholpseudoinv4x4f(void* pScr,
    complex_float *x,
    const complex_float * A,
    const float32_t  sigma2)
{
    complex_float * restrict D; int SD;
    complex_float * restrict R; int SR;
    complex_float * restrict y; int SY;

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
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);
    R = (complex_float *)pScr;
    D = (complex_float *)(((uintptr_t)R) + SR);
    y = (complex_float *)(((uintptr_t)D) + SD);
    pScr = (complex_float *)(((uintptr_t)y) + SY);

    cmatcholdecomp4x4f(pScr, R, D, A, sigma2);
    fwd4x4f(y, R, D, A);
    bkw4x4f(x, R, D, y);
}

size_t  cmatcholpseudoinv4x4f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY/*, SRt*/;

    SD = 4;
    SR = (((4 + 1) * 4) >> 1);
    SY = 4 * 4;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);

    s_dc = cmatcholdecomp4x4f_getScratchSize();
    return SD + SR + SY + s_dc;
}
#else
DISCARD_FUN(void, cmatcholpseudoinv4x4f, (void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t sigma2))

size_t  cmatcholpseudoinv4x4f_getScratchSize()
{
	return 0;
}
#endif
