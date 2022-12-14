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

#if (HAVE_VFPU || HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)

/* Load 32bit and replicate to xtfloatx2*/
#define _L32_SX2_IP(a,b,c) \
{ \
ae_int32x2 tmp; \
AE_L32_IP(tmp, castxcc(ae_int32, b), c); \
a = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp); }

/*-------------------------------------------------------------------------
Cholesky decomposition for Pseudo-Inversion
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem:
A[MxN]*X[NxP]=B[MxP], MxN.
The decomposition results in an upper triangular complex matrix R[NxN]
with real and positive numbers on the main diagonal, such that
R'*R = A'*A + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]          matrix A. For fixed point API, the representation is Q(qA)
sigma2          regularization term. For fixed point, the representation 
                should be Q(2*qA-30)
qRA             qR-qA; difference between fixed point representations of R
                and A matrices (for the fixed point API only). Should be 
                equal or less than 0 (typically -2).
Output:
R[((N+1)*N)/2]  upper triangular matrix R
D[N]            reciprocals of the main diagonal. NOTE: for the fixed point API,
                these data are stored internally in special format with separate
                mantissa and exponent for better accuracy and dynamic range 
                control. So, even for the real data, they stored as pairs of 2
                integers and packed to the complex_fract32 format
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
8-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholdecomp8x8f(void * pScr,
    complex_float * R,
    complex_float * D,
    const complex_float * A,
    const float32_t       sigma2)
{
    static const f_cplxcholFwdrecf f_Fwdrec[8] = { cplxcholFwdrec1f, cplxcholFwdrec2f, cplxcholFwdrec3f, cplxcholFwdrec4f, cplxcholFwdrec5f, cplxcholFwdrec6f, cplxcholFwdrec7f, cplxcholFwdrec8f };
    int n;
    complex_float* y = R; // pointer to the new column in R
    complex_float* prep = (complex_float*)pScr;
    NASSERT(R);
    NASSERT(D);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    // compute A'*A
    cmatcholpreprocess8x8f(NULL, prep, A, sigma2);

    for (n = 0; n<8; n++)
    {
        y += n; // go to the next column
        // make forward recursion to update n new column elements
        f_Fwdrec[n](y, R, D, prep + n, 8);
        // update n-th diagonal element
        cplxcholDiagUpdf(y, D + n, prep + n*9, n);
    }
}

size_t  cmatcholdecomp8x8f_getScratchSize()
{
    return 8 * 8 * SZ_CF32;
}


#else
DISCARD_FUN(void, cmatcholdecomp8x8f, (void * pScr,
    complex_float * R,
    complex_float * D,
    const complex_float * A,
    const float32_t sigma2))

size_t  cmatcholdecomp8x8f_getScratchSize()
{
    return 0;
}
#endif
