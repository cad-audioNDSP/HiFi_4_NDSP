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
#include "chol32x32_common.h"

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
void  cmatcholpseudoinv6x6_32x32(void* pScr,
        complex_fract32 *x,
        const complex_fract32 * A,
        const int32_t  sigma2, 
        int qRA, int qYBRA, int qXYR)
{
    const int N=6;
    tcmatcholpseudoinv_32x32_Scratch scr;
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);

    tcmatcholpseudoinv_32x32_getScratch(&scr, pScr, N);
    pScr=scr.pEnd;
    cmatcholdecomp6x6_32x32(pScr,scr.R,scr.D,A,sigma2,qRA);
    cmatcholtransformR(scr.Rt,scr.R,N);
    cmatcholfwdsubstnxn(scr.y,scr.R,scr.D,A,N,qYBRA);
    cmatcholbkwsubstnxn(x ,scr.Rt,scr.D,scr.y, N, qXYR);
}

size_t  cmatcholpseudoinv6x6_32x32_getScratchSize()
{ 
    size_t sz;
    tcmatcholpseudoinv_32x32_Scratch scr;
    sz= cmatcholdecomp6x6_32x32_getScratchSize();
    tcmatcholpseudoinv_32x32_getScratch(&scr, NULL, 6);
    return sz+scr.sz;
}
