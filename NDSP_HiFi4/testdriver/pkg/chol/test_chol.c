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
/*
 * Test procedures for Cholesky decomposition (block data)
 */

#include "NatureDSP_types.h"
#include "NatureDSP_Math.h"
#include "NatureDSP_Signal_matinv.h"
#include "malloc16.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
/* Test Engine API. */
#include "testeng.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define LOG2_SIMD_WIDTH 2

// get allocated space per one matrix
// for one matrix we won't insert holes so allocated space equals number of elements in matrix
static int getSpace(int S,int szElem)
{
    return S;   /* PDX API: no holes between matrices in single-matrix API!!! */
}

/* chol API flags */
#define CHOLAPI_FIX_SZ        0x0001 /* flag for fixed-size API like chol8x8n            */
#define CHOLAPI_FLT_SZ        0x0002 /* flag for floating-size API like cholmxnn         */
#define CHOLAPI_SZTYPE_MASK   0x0003

#define CHOLAPI_PACKED        0x0004 /* flag for block-ordered functions like cholmxnn  */
#define CHOLAPI_STREAM        0x0008 /* flag for stream-ordered functions like cholmxns */
#define CHOLAPI_LAYOUT_MASK   0x000C
#define CHOLAPI_USESIGMA2     0x0010 /* flag for functions that use sigma2 regularization term */

static void cholmxnn_32b_unpackD(void * dst, const void * src, int N, int L)
{
    int32_t *D = (int32_t *)dst;
    const float64_t *X = (const float64_t *)src;
    int l,n,Nd;
    Nd = getSpace(N, sizeof(complex_fract32));
    for (l=0; l<L; l++)
    {
        for (n=0; n<N; n++)
        {
            float64_t mant;
            int exp;

            mant = frexp(X[l*Nd+n],&exp);/* q[-exp] */
            /* treat D as Q[62 - exp] */
            D[2*(l*Nd+n)+1] = (int32_t)exp+31;
            D[2*(l*Nd+n)+0] = (int32_t)MIN(MAX_INT32,round(ldexp(mant, 31)));
        }
    }
}

static void cholmxnn_32b_packD(void * dst, const void * src, int N, int L)
{
    float64_t *X = (float64_t *)dst;
    const int32_t *D = (const int32_t *)src;
    int l,n,Nd;
    Nd = getSpace(N, sizeof(complex_fract32));
    for (l=0; l<L; l++)
    {
        for (n=0; n<Nd; n++)
        {
            float64_t mant;
            int exp;
            mant= D[2*(l*Nd+n)+0];/* q[62 - exp] */
            exp = D[2*(l*Nd+n)+1];
            X[l*Nd+n]=ldexp(mant, exp - 62);
        }
    }
}

/* API definition structure. */
typedef struct 
{
  size_t (*cholScratchSz)(); 
  size_t (*fwdScratchSz)(); 
  size_t (*bkwScratchSz)();        
  size_t (*mmseScratchSz)(); 
  size_t (*preprocessScratchSz)(); 
  size_t (*pinvScratchSz)(); 

    /* Helping function, convert reciprocal of diagonal elements from the format used in SEQ-file into the format used in function */
    void (*te_cholFxn_unpackD)(void * dst, const void * src, int N, int L);
    /* Helping function, convert reciprocal of diagonal elements from the format used in function into the format used in SEQ-file */
    void (*te_cholFxn_packD)(void * dst, const void * src, int N, int L);

} tCholApi;

static const tCholApi cmatchol4x4_32x32_Api   =    { cmatcholdecomp4x4_32x32_getScratchSize,   cmatcholfwdsubst4x4_32x32_getScratchSize  ,cmatcholbkwsubst4x4_32x32_getScratchSize   ,cmatcholmmsesolver4x4_32x32_getScratchSize   ,cmatcholpreprocess4x4_32x32_getScratchSize  ,cmatcholpseudoinv4x4_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi cmatchol6x6_32x32_Api   =    { cmatcholdecomp6x6_32x32_getScratchSize,   cmatcholfwdsubst6x6_32x32_getScratchSize  ,cmatcholbkwsubst6x6_32x32_getScratchSize   ,cmatcholmmsesolver6x6_32x32_getScratchSize   ,cmatcholpreprocess6x6_32x32_getScratchSize  ,cmatcholpseudoinv6x6_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi cmatchol8x8_32x32_Api   =    { cmatcholdecomp8x8_32x32_getScratchSize,   cmatcholfwdsubst8x8_32x32_getScratchSize  ,cmatcholbkwsubst8x8_32x32_getScratchSize   ,cmatcholmmsesolver8x8_32x32_getScratchSize   ,cmatcholpreprocess8x8_32x32_getScratchSize  ,cmatcholpseudoinv8x8_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi cmatchol10x10_32x32_Api =    { cmatcholdecomp10x10_32x32_getScratchSize, cmatcholfwdsubst10x10_32x32_getScratchSize,cmatcholbkwsubst10x10_32x32_getScratchSize ,cmatcholmmsesolver10x10_32x32_getScratchSize ,cmatcholpreprocess10x10_32x32_getScratchSize,cmatcholpseudoinv10x10_32x32_getScratchSize, cholmxnn_32b_unpackD, cholmxnn_32b_packD };

static const tCholApi  matchol4x4_32x32_Api   =    {  matcholdecomp4x4_32x32_getScratchSize,   matcholfwdsubst4x4_32x32_getScratchSize   ,matcholbkwsubst4x4_32x32_getScratchSize    ,matcholmmsesolver4x4_32x32_getScratchSize    ,cmatcholpreprocess4x4_32x32_getScratchSize  , matcholpseudoinv4x4_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi  matchol6x6_32x32_Api   =    {  matcholdecomp6x6_32x32_getScratchSize,   matcholfwdsubst6x6_32x32_getScratchSize   ,matcholbkwsubst6x6_32x32_getScratchSize    ,matcholmmsesolver6x6_32x32_getScratchSize    ,cmatcholpreprocess6x6_32x32_getScratchSize  , matcholpseudoinv6x6_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi  matchol8x8_32x32_Api   =    {  matcholdecomp8x8_32x32_getScratchSize,   matcholfwdsubst8x8_32x32_getScratchSize   ,matcholbkwsubst8x8_32x32_getScratchSize    ,matcholmmsesolver8x8_32x32_getScratchSize    ,cmatcholpreprocess8x8_32x32_getScratchSize  , matcholpseudoinv8x8_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi  matchol10x10_32x32_Api =    {  matcholdecomp10x10_32x32_getScratchSize, matcholfwdsubst10x10_32x32_getScratchSize ,matcholbkwsubst10x10_32x32_getScratchSize  ,matcholmmsesolver10x10_32x32_getScratchSize  ,cmatcholpreprocess10x10_32x32_getScratchSize, matcholpseudoinv10x10_32x32_getScratchSize, cholmxnn_32b_unpackD, cholmxnn_32b_packD };

static const tCholApi chol4x4x1nfApi =    { cmatcholdecomp4x4f_getScratchSize,   cmatcholfwdsubst4x4f_getScratchSize,   cmatcholbkwsubst4x4f_getScratchSize,   cmatcholmmsesolver4x4f_getScratchSize,   cmatcholpreprocess4x4f_getScratchSize,   cmatcholpseudoinv4x4f_getScratchSize };
static const tCholApi chol6x6x1nfApi =    { cmatcholdecomp6x6f_getScratchSize,   cmatcholfwdsubst6x6f_getScratchSize,   cmatcholbkwsubst6x6f_getScratchSize,   cmatcholmmsesolver6x6f_getScratchSize,   cmatcholpreprocess6x6f_getScratchSize,   cmatcholpseudoinv6x6f_getScratchSize };
static const tCholApi chol8x8x1nfApi =    { cmatcholdecomp8x8f_getScratchSize,   cmatcholfwdsubst8x8f_getScratchSize,   cmatcholbkwsubst8x8f_getScratchSize,   cmatcholmmsesolver8x8f_getScratchSize,   cmatcholpreprocess8x8f_getScratchSize,   cmatcholpseudoinv8x8f_getScratchSize };
static const tCholApi chol10x10x1nfApi =  { cmatcholdecomp10x10f_getScratchSize, cmatcholfwdsubst10x10f_getScratchSize, cmatcholbkwsubst10x10f_getScratchSize, cmatcholmmsesolver10x10f_getScratchSize, cmatcholpreprocess10x10f_getScratchSize, cmatcholpseudoinv10x10f_getScratchSize };
static const tCholApi rchol4x4x1nfApi =   { matcholdecomp4x4f_getScratchSize,    matcholfwdsubst4x4f_getScratchSize,    matcholbkwsubst4x4f_getScratchSize,    matcholmmsesolver4x4f_getScratchSize,    matcholpreprocess4x4f_getScratchSize,    matcholpseudoinv4x4f_getScratchSize };
static const tCholApi rchol6x6x1nfApi =   { matcholdecomp6x6f_getScratchSize,    matcholfwdsubst6x6f_getScratchSize,    matcholbkwsubst6x6f_getScratchSize,    matcholmmsesolver6x6f_getScratchSize,    matcholpreprocess6x6f_getScratchSize,    matcholpseudoinv6x6f_getScratchSize };
static const tCholApi rchol8x8x1nfApi =   { matcholdecomp8x8f_getScratchSize,    matcholfwdsubst8x8f_getScratchSize,    matcholbkwsubst8x8f_getScratchSize,    matcholmmsesolver8x8f_getScratchSize,    matcholpreprocess8x8f_getScratchSize,    matcholpseudoinv8x8f_getScratchSize };
static const tCholApi rchol10x10x1nfApi = { matcholdecomp10x10f_getScratchSize,  matcholfwdsubst10x10f_getScratchSize,  matcholbkwsubst10x10f_getScratchSize,  matcholmmsesolver10x10f_getScratchSize,  matcholpreprocess10x10f_getScratchSize,  matcholpseudoinv10x10f_getScratchSize };
/*
static const tCholApi cholmxnxpnfApi = { cmatcholf_getScratchSize, cmatcholfwdf_getScratchSize, cmatcholbkwf_getScratchSize, cmatcholmmsef_getScratchSize };
static const tCholApi rcholmxnxpnfApi = { matcholf_getScratchSize, matcholfwdf_getScratchSize, matcholbkwf_getScratchSize, matcholmmsef_getScratchSize };
static const tCholApi cholmxnxpnApi  = {cmatchol_getScratchSize  ,cmatcholfwd_getScratchSize  ,cmatcholbkw_getScratchSize ,cmatcholmmse_getScratchSize };
static const tCholApi rcholmxnxpnApi  = {matchol_getScratchSize  ,matcholfwd_getScratchSize  ,matcholbkw_getScratchSize ,matcholmmse_getScratchSize };
*/

typedef struct
{
    int NA;
    int NR;
    int ND;
    int NY;
    int NB;
    int NX;
    int NQ;
    int NSigma2;
} tCholNelems;

/* Compute number of elements needed to allocate for Cholesky decomposition functions */
static tCholNelems te_szHelper_chol( tTestEngContext * context )
{
    int M,N/*,L*/;
    tCholNelems Nelem;

    ASSERT( context && context->seqFile && context->desc->extraPtr );

    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
//    L = MAX(0,context->args.dim[2]);

    memset( &Nelem, 0, sizeof(Nelem) );
    /* define the number of elements for input/output data */
    if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_PACKED)
    {
        int sz=0, szD=0;

        switch(context->desc->fmt & FMT_DTYPE_MASK)
        {
        case FMT_FRACT16:
            sz = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract16) : sizeof(int16_t);
            szD = sizeof(complex_fract16);
            break;
        case FMT_FRACT32:
            sz = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract32) : sizeof(int32_t);
            szD = sizeof(complex_fract32);
            break;
        case FMT_FLOAT32:
            sz = szD = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_float) : sizeof(float32_t);
            break;
        default: ASSERT(0);
        }
//        if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT16)
//        {
//            sz = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract16) : sizeof(int16_t);
//            szD = sizeof(complex_fract16);
//        }
//        else if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT32)
//        {
//            sz = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract32) : sizeof(int32_t);
//            szD = sizeof(complex_fract32);
//        }
//        else
//        {
//            sz = szD = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_float) : sizeof(float32_t);
//        }
        Nelem.NR = getSpace((N*(N+1))/2, sz);
        Nelem.ND = getSpace(N, szD);
        Nelem.NA = getSpace(M*N, sz);
    }
    else if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_STREAM)
    {
        Nelem.NR = N*N;
        Nelem.NA = M*N;
        Nelem.ND = N;
    }
    else
    {
        ASSERT(0);
    }

    /* define if used term sigma2 or not */
    if (context->desc->extraParam & CHOLAPI_USESIGMA2)
    {
        Nelem.NSigma2 = 1;
    }
    else
    {
        Nelem.NSigma2 = 0;
    }

    /* define the number of fract parameters Q */
    if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT32 || (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT16)/* all floating-point funcs */
    {
        Nelem.NQ = 0;
    }
    else if (context->desc->fmt == (FMT_CPLX|FMT_FRACT16))/* chols, choln */
    {
        Nelem.NQ = 4;/* qA, qB, qX, qY*/
    }
    else/* rchols, dchols, dcholn, chols_32x32, choln_32x32 */
    {
        Nelem.NQ = 5;/* qA, qB, qR, qY, qX*/
    }

    return Nelem;
} /* te_szHelper_chol() */



/* Allocate vectors and load the data set for Cholesky decomposition functions */
int te_loadFxn_chol( tTestEngContext * context )
{
    int /*M,N,*/L,Nr,Nd,Nq,Na,NSigma2,Narg,NTotalArgs,res = 0;
    int fmtD=0, fmtSigma2=0;
    tCholNelems Nelems;

    ASSERT( context && context->seqFile && context->desc->extraPtr );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );

    /* Get input data size and format. */
//    M = MAX(0,context->args.dim[0]);
//    N = MAX(0,context->args.dim[1]);
    L = MAX(0,context->args.dim[2]);
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
    case FMT_FRACT32:
        fmtD     =FMT_REAL|FMT_FLOAT64;
        fmtSigma2=FMT_REAL|FMT_FRACT32;
        break;
    case FMT_FRACT16:
        fmtD     =FMT_REAL|FMT_FLOAT32;
        fmtSigma2=FMT_REAL|FMT_INT32;
        break;
    case FMT_FLOAT32:
        fmtD     =FMT_REAL|FMT_FLOAT32;
        fmtSigma2=FMT_REAL|FMT_FLOAT32;
        break;
    case FMT_FLOAT16:
        fmtD     =FMT_REAL|FMT_FLOAT16;
        fmtSigma2=FMT_REAL|FMT_FLOAT16;
        break;
    }
//    fmtD = (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT32 ? FMT_REAL|FMT_FLOAT64 : FMT_REAL|FMT_FLOAT32;
//    fmtSigma2 = (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT16 ? FMT_REAL|FMT_INT32 : FMT_REAL|FMT_FLOAT32;

    Nelems = te_szHelper_chol(context);
    Nr=Nelems.NR;
    Na=Nelems.NA;
    Nd=Nelems.ND;
    Nq=Nelems.NQ;
    NSigma2 = Nelems.NSigma2;

    /* Allocate input data vector memory. */
    NTotalArgs = 0;
    Narg =vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X  , Na*L, /* A */
                    &context->dataSet.Z  , Nr*L, /* R */
                    &context->dataSet.Zlo, Nr*L, /* R */
                    &context->dataSet.Zhi, Nr*L, /* R */
                    &context->dataSet.auxVec[1], 2*Nd*L, /* D from output */
                    0);
    NTotalArgs += 5;
    Narg+=vecsAlloc( context->desc->isAligned, fmtD,
                    &context->dataSet.W  , Nd*L, /* D from SEQ-file */
                    &context->dataSet.Wlo, Nd*L, /* D from SEQ-file */
                    &context->dataSet.Whi, Nd*L, /* D from SEQ-file */
                    0 );
    NTotalArgs += 3;

    if (Nq > 0)/* allocate memory for parameters Q */
    {
        Narg+=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.Y  , Nq, /* Q */
                        0);
        NTotalArgs += 1;
    }
    if (NSigma2 > 0)/* allocate memory for sigma2 */
    {
        Narg+=vecsAlloc( context->desc->isAligned, fmtSigma2,
                        &context->dataSet.U  , NSigma2*L, /* sigma2 */
                        0);
        NTotalArgs += 1;
    }

    if(Narg!=NTotalArgs)
    {
        printf( "te_loadFxn_chol(): failed to allocate memory\n" );
    }
    else
    {
        /* clean output of function just to protect from ferret warnings appearing 
            when test environment reads data in the holes between matrices */
        memset(vecGetElem(&context->dataSet.Z,0),0,vecGetSize(&context->dataSet.Z));
        memset(vecGetElem(&context->dataSet.W,0),0,vecGetSize(&context->dataSet.W));
        memset(vecGetElem(&context->dataSet.auxVec[1],0),0,vecGetSize(&context->dataSet.auxVec[1]));

        /* Load vectors data from the SEQ-file. */
        res = 1;
        /* sigma2 - for cholnf, rcholnf */
//        if ((NSigma2 > 0) && ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_PACKED))
//        {
//            res &= seqFileReadVec(context->seqFile, &context->dataSet.U);
//        }
        /* A */
        res &= seqFileReadVec(context->seqFile, &context->dataSet.X);
        /* sigma2 */
        if ((NSigma2 > 0) )
        {
            res &= seqFileReadVec(context->seqFile, &context->dataSet.U);
        }
        /* Q */
        if (Nq > 0)
        {
            res &= seqFileReadVec(context->seqFile, &context->dataSet.Y);
        }
        /* Rlo, Rhi, Dlo, Dhi */
        res &= seqFileReadVecs(context->seqFile,
                               &context->dataSet.Zlo,
                               &context->dataSet.Zhi,
                               &context->dataSet.Wlo,
                               &context->dataSet.Whi,
                               0);

        if ( !res )
        {
            printf( "te_loadFxn_chol(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_chol() */

/* processing function for Cholesky decomposition functions */
void te_processFxn_chol( tTestEngContext * context )
{
    typedef void tFnProcess_sigma2_blk_int32(void *pScr, void *R, void *D, const void *A, int32_t sigma2, int qRA);
    int32_t N, M, L;
    void *A,*R,*D;
    void *w;// reciprocal of main diagonal in floating point representation
    const tCholApi *pCholContext;
    te_fun_ptr_t fxn;
    fxn = context->target.fut;

    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    L = MAX(0,context->args.dim[2]);
    (void)M,(void)N,(void)L;
    A = vecGetElem(&context->dataSet.X,0);
    R = vecGetElem(&context->dataSet.Z,0);
    w = vecGetElem(&context->dataSet.W,0);
    D = vecGetElem(&context->dataSet.auxVec[1],0);
    pCholContext = (const tCholApi *)context->desc->extraPtr;

    /* Apply function */
    tVec vScr;
    int scr_sz;

    /* Allocate scratch buffer */
    scr_sz = pCholContext->cholScratchSz();
    if (!vecAlloc(&vScr, scr_sz, context->desc->isAligned, FMT_REAL|FMT_UINT8, NULL)) 
    {
        printf( "te_processFxn_chol(): failed to allocate memory\n");
    }
    else
    {
        int qA, qR, qRA;
        void *pScr;
        void *sigma2;

        pScr = vecGetElem(&vScr, 0);
        {   /* reporting */
            tReportFUT fut[2];
            fut[0]=(tReportFUT)context->target.fut;
            fut[1]=(tReportFUT)pCholContext->cholScratchSz;
            vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
        }
        switch(context->desc->extraParam)
        {
            case CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2:
                qA = *vecGetElem_i32(&context->dataSet.Y, 0);
                qR = *vecGetElem_i32(&context->dataSet.Y, 2);
                qRA = qR - qA;
                sigma2 = vecGetElem(&context->dataSet.U,0);
                ((tFnProcess_sigma2_blk_int32*)fxn)(pScr,R,D,A,((int32_t*)sigma2)[0],qRA);
                break;
            default:
                ASSERT(0);
        }
        vecFree(&vScr);
    }

    /* copy D to w with right formatting */
    pCholContext->te_cholFxn_packD(w, D, N, 1);
} /* te_processFxn_chol() */


/* Compute number of elements needed to allocate for Cholesky backward recursion functions */
static tCholNelems te_szHelper_cholbkw( tTestEngContext * context )
{
    int N,P/*,L*/;
    tCholNelems Nelem;

    ASSERT( context && context->seqFile && context->desc->extraPtr );

    N = MAX(0,context->args.dim[0]);
    P = MAX(0,context->args.dim[1]);
//    L = MAX(0,context->args.dim[2]);

    memset( &Nelem, 0, sizeof(Nelem) );

    /* define the number of elements for input/output data */
    if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_PACKED)
    {
        int sz=0, szD=0;
        switch(context->desc->fmt & FMT_DTYPE_MASK)
        {
        case FMT_FRACT16:
            sz  = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract16) : sizeof(int16_t);
            szD = sizeof(complex_fract16);
            break;
        case FMT_FRACT32:
            sz = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract32) : sizeof(int32_t);
            szD = sizeof(complex_fract32);
            break;
        case FMT_FLOAT32:
            sz = szD = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_float) : sizeof(float32_t);
            break;
        default: ASSERT(0);
        }
//        if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT16)
//        {
//            sz  = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract16) : sizeof(int16_t);
//            szD = sizeof(complex_fract16);
//        }
//        else if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT32)
//        {
//            sz = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract32) : sizeof(int32_t);
//            szD = sizeof(complex_fract32);
//        }
//        else
//        {
//            sz = szD = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_float) : sizeof(float32_t);
//        }

        Nelem.NR = getSpace((N*(N + 1))/2, sz);
        Nelem.ND = getSpace(N, szD);
        Nelem.NY = getSpace(N*P, sz);
        Nelem.NX = getSpace(N*P, sz);
    }
    else if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_STREAM)
    {
        Nelem.NR = N*N;
        Nelem.ND = N;
        Nelem.NX = N*P;
        Nelem.NY = N*P;
    }
    else
    {
        ASSERT(0);
    }

    /* define the number of fract parameters Q */
    if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT32 || (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT16)
    {
        Nelem.NQ = 0;
    }
    else if (context->desc->fmt == (FMT_CPLX|FMT_FRACT16))/* chols, choln */
    {
        Nelem.NQ = 4;/* qA, qB, qX, qY*/
    }
    else
    {
        Nelem.NQ = 5;/* qA, qB, qR, qY, qX*/
    }

    return Nelem;
} /* te_szHelper_cholbkw() */

/* Allocate vectors and load the data set for Cholesky backward recursion functions */
int te_loadFxn_cholbkw( tTestEngContext * context )
{
    int /*N,P,*/L,res = 0;
    int Nx,Ny,Nr,Nd,Nq,Narg,NTotalArgs;
    int fmtD=0;
    tCholNelems Nelems;

    ASSERT( context && context->seqFile && context->desc->extraPtr );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    /* Get input data size and format. */
//    N = MAX(0,context->args.dim[0]);
//    P = MAX(0,context->args.dim[1]);
    L = MAX(0,context->args.dim[2]);

    Nelems = te_szHelper_cholbkw(context);
    Nr = Nelems.NR;
    Nd = Nelems.ND;
    Ny = Nelems.NY;
    Nx = Nelems.NX;
    Nq = Nelems.NQ;
//    fmtD = (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT32 ? FMT_REAL|FMT_FLOAT64 : FMT_REAL|FMT_FLOAT32;
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
    case FMT_FRACT16:
        fmtD = FMT_REAL|FMT_FLOAT32;
        break;
    case FMT_FRACT32:
        fmtD = FMT_REAL|FMT_FLOAT64;
        break;
    case FMT_FLOAT16:
        fmtD =FMT_REAL|FMT_FLOAT16;
        break;
    case FMT_FLOAT32:
        fmtD =FMT_REAL|FMT_FLOAT32;
        break;
    default: ASSERT(0);
    }

    /* Allocate input data vector memory. */
    Narg =vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X,         Nr*L, // R
                    &context->dataSet.U,         Ny*L, // Y
                    &context->dataSet.Z,         Nx*L, // X
                    &context->dataSet.Zlo,       Nx*L, // X
                    &context->dataSet.Zhi,       Nx*L, // X
                    &context->dataSet.auxVec[1], 2*Nd*L, // D to input
                    0);
    NTotalArgs = 6;
    Narg+=vecsAlloc( context->desc->isAligned, fmtD,
                    &context->dataSet.Y  , Nd*L,     // D from SEQ-file
                    0);
    NTotalArgs += 1;
    if (Nq > 0)
    {
        Narg+=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.auxVec[0]  , Nq,
                        0);
        NTotalArgs += 1;
    }

    if(Narg!=NTotalArgs)
    {
        printf( "te_loadFxn_cholbkw(): failed to allocate memory\n" );
    }
    else
    {
        /* clean output of function just to protect from ferret warnings appearing 
           when test environment reads data in the holes between matrices */
        memset(vecGetElem(&context->dataSet.Z,0),0,vecGetSize(&context->dataSet.Z));
        memset(vecGetElem(&context->dataSet.auxVec[1],0),0,vecGetSize(&context->dataSet.auxVec[1]));

        /* Load vectors data from the SEQ-file. */
        res = 1;
        if (Nq > 0)
        {
            res &= seqFileReadVec(context->seqFile, &context->dataSet.auxVec[0]);
        }
        res &= seqFileReadVecs(context->seqFile,
                               &context->dataSet.X,
                               &context->dataSet.Y,
                               &context->dataSet.U,
                               &context->dataSet.Zlo,
                               &context->dataSet.Zhi,
                               0);

        if (!res)
        {
            printf( "te_loadFxn_cholbkw(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholbkw() */

/* processing function for Cholesky backward recursion functions */
void te_processFxn_cholbkw( tTestEngContext * context )
{
    typedef void tFnProcess(void *pScr, void *x, const void *R, const void *D, const void *y, int qXYR);
    typedef void tFnProcessf_blk(void *pScr, void *x, const void *R, const void *D, const void *y);
    int32_t N,P,L,Nq;
    int qR,qX,qY,qXYR;
    void *R, *D, *X, *Y;
    void *w;// reciprocal of main diagonal in floating point representation
    const tCholApi *pCholContext;
    te_fun_ptr_t fxn;
    fxn = context->target.fut;

    N = MAX(0,context->args.dim[0]);
    P = MAX(0,context->args.dim[1]);
    L = MAX(0,context->args.dim[2]);
    (void)N,(void)P,(void)L;
    R = vecGetElem(&context->dataSet.X, 0);
    Y = vecGetElem(&context->dataSet.U, 0);
    X = vecGetElem(&context->dataSet.Z, 0);
    D = vecGetElem(&context->dataSet.auxVec[1], 0);
    w = vecGetElem(&context->dataSet.Y,0);

    qXYR = 0;
    Nq = te_szHelper_cholbkw(context).NQ;
    if (Nq > 0)
    {
        int32_t *Q;
        Q = ((int32_t *)vecGetElem(&context->dataSet.auxVec[0], 0));
        if (Nq == 4)
        {
            /* qA, qB, qX, qY */
            qR = Q[0];
            qX = Q[2];
            qY = Q[3];
            qXYR = qX - qY + qR;
        }
        else
        {
            /* qA, qB, qR, qY, qX */
            qR = Q[2];
            qY = Q[3];
            qX = Q[4];
            qXYR = qX - qY + qR;
        }
    }
    pCholContext = (const tCholApi *)context->desc->extraPtr;
    /* Copy w to D with right formatting */
    pCholContext->te_cholFxn_unpackD(D, w, N, L);

    /* Apply function */
    {
        tVec vScr;
        int scr_sz;

        /* Allocate scratch buffer */
        scr_sz = pCholContext->bkwScratchSz();
        if (!vecAlloc(&vScr, scr_sz, context->desc->isAligned, FMT_REAL|FMT_UINT8, NULL)) 
        {
            printf( "te_processFxn_cholbkw(): failed to allocate memory\n");
        }
        else
        {
            void * pScr;
            pScr = vecGetElem(&vScr, 0);
            {    /* reporting */
                tReportFUT fut[2];
                fut[0]=(tReportFUT)context->target.fut;
                fut[1]=(tReportFUT)pCholContext->bkwScratchSz;
                vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
            }
            switch(context->desc->extraParam)
            {
                case CHOLAPI_FIX_SZ|CHOLAPI_PACKED:
                    if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT32 || (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT16)
                        ((tFnProcessf_blk*)fxn)(pScr, X, R, D, Y);
                    else
                        ((tFnProcess*)fxn)(pScr, X, R, D, Y, qXYR);
                    break;
                default: ASSERT(0);
            }
            vecFree(&vScr);
        }
    }
} /* te_processFxn_cholbkw() */

static tCholNelems te_szHelper_cholfwd( tTestEngContext * context )
{
    int M,N/*,L*/,P;
    tCholNelems Nelem;

    ASSERT( context && context->seqFile && context->desc->extraPtr );

    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    P = MAX(0,context->args.dim[2]);
//    L = MAX(0,context->args.dim[3]);

    memset( &Nelem, 0, sizeof(Nelem) );
    if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_PACKED)
    {
        int sz=0, szD=0;

        switch(context->desc->fmt & FMT_DTYPE_MASK)
        {
        case FMT_FRACT16:
            sz =(context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract16) : sizeof(int16_t);;
            szD=sizeof(complex_fract16);
            break;
        case FMT_FRACT32:
            sz =(context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract32) : sizeof(int32_t);;
            szD=sizeof(complex_fract32);
            break;
        case FMT_FLOAT32:
            sz = szD = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_float) : sizeof(float32_t);
            break;
        default: ASSERT(0);
        }
//        if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT16)
//        {
//            sz  = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract16) : sizeof(int16_t);
//            szD = sizeof(complex_fract16);
//        }
//        else if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT32)
//        {
//            sz  = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_fract32) : sizeof(int32_t);
//            szD = sizeof(complex_fract32);
//        }
//        else
//        {
//            sz = szD = (context->desc->fmt & FMT_CPLX) ? sizeof(complex_float) : sizeof(float32_t);
//        }

        Nelem.NR = getSpace((N*(N+1))/2, sz);
        Nelem.ND = getSpace(N, szD);
        Nelem.NA = getSpace(M*N, sz);
        Nelem.NY = getSpace(N*P, sz);
        Nelem.NB = getSpace(M*P, sz);
    }
    else if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_STREAM)
    {
        Nelem.NR = N*N;
        Nelem.ND = N;
        Nelem.NA = M*N;
        Nelem.NY = N*P;
        Nelem.NB = M*P;
    }
    else
    {
        ASSERT(0);
    }

    /* define the number of fract parameters Q */
    if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT32 || (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT16)
    {
        Nelem.NQ = 0;
    }
    else if (context->desc->fmt == (FMT_CPLX|FMT_FRACT16))/* chols, choln */
    {
        Nelem.NQ = 4;/* qA, qB, qX, qY*/
    }
    else
    {
        Nelem.NQ = 5;/* qA, qB, qR, qY, qX*/
    }

    return Nelem;
} /* te_szHelper_cholfwd() */
/* Allocate vectors and load the data set for Cholesky forward substitution functions */
int te_loadFxn_cholfwd( tTestEngContext * context )
{
    int /*M,N,P,*/L,res = 0;
    int Nr,Nd,Nq,Na,Nb,Ny,Narg,NTotalArgs;
    int fmtD=0;
    tCholNelems Nelems;
    ASSERT( context && context->seqFile && context->desc->extraPtr );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );

    /* Get input data size and format. */
//    M = MAX(0,context->args.dim[0]);
//    N = MAX(0,context->args.dim[1]);
//    P = MAX(0,context->args.dim[2]);
    L = MAX(0,context->args.dim[3]);

    Nelems = te_szHelper_cholfwd(context);
    Nr = Nelems.NR;
    Nd = Nelems.ND;
    Ny = Nelems.NY;
    Na = Nelems.NA;
    Nb = Nelems.NB;
    Nq = Nelems.NQ;
//    fmtD = (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FRACT32 ? FMT_REAL|FMT_FLOAT64 : FMT_REAL|FMT_FLOAT32;
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
    case FMT_FRACT16:
        fmtD = FMT_REAL|FMT_FLOAT32;
        break;
    case FMT_FRACT32:
        fmtD = FMT_REAL|FMT_FLOAT64;
        break;
    case FMT_FLOAT16:
        fmtD =FMT_REAL|FMT_FLOAT16;
        break;
    case FMT_FLOAT32:
        fmtD =FMT_REAL|FMT_FLOAT32;
        break;
    default: ASSERT(0);
    }

    /* Allocate input data vector memory. */
    Narg =vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X,         Nr*L, // R
                    &context->dataSet.U,         Na*L, // A
                    &context->dataSet.V,         Nb*L, // B
                    &context->dataSet.Z,         Ny*L, // Y
                    &context->dataSet.Zlo,       Ny*L, // Y
                    &context->dataSet.Zhi,       Ny*L, // Y
                    &context->dataSet.auxVec[1], 2*Nd*L, // D to input
                    0);
    NTotalArgs = 7;
    Narg+=vecsAlloc( context->desc->isAligned, fmtD,
                    &context->dataSet.Y  , Nd*L,        // D from SEQ-file
                    0);
    NTotalArgs += 1;
    if (Nq > 0)
    {
        Narg+=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.auxVec[0], Nq,
                        0);
        NTotalArgs += 1;
    }
    if(Narg!=NTotalArgs)
    {
        printf( "te_loadFxn_cholfwd(): failed to allocate memory\n" );
    }
    else
    {
        /* clean output of function just to protect from ferret warnings appearing 
           when test environment reads data in the holes between matrices */
        memset(vecGetElem(&context->dataSet.Z,0),0,vecGetSize(&context->dataSet.Z));
        memset(vecGetElem(&context->dataSet.auxVec[1],0),0,vecGetSize(&context->dataSet.auxVec[1]));

        /* Load vectors data from the SEQ-file. */
        if (Nq > 0)
        {
            res = seqFileReadVecs(context->seqFile,
                &context->dataSet.auxVec[0],
                &context->dataSet.X,
                &context->dataSet.Y,
                &context->dataSet.U,
                &context->dataSet.V,
                &context->dataSet.Zlo,
                &context->dataSet.Zhi,
                0);
        }
        else
        {
            res = seqFileReadVecs(context->seqFile,
                &context->dataSet.X,
                &context->dataSet.Y,
                &context->dataSet.U,
                &context->dataSet.V,
                &context->dataSet.Zlo,
                &context->dataSet.Zhi,
                0);
        }

        if (!res)
        {
            printf( "te_loadFxn_cholfwd(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholfwd() */
/* processing function for Cholesky forward substitution functions */
void te_processFxn_cholfwd( tTestEngContext * context )
{
    typedef void tFnProcess(void *pScr, void *y, const void *R, const void *D, const void *A, const void *B, int qYBRA);
    typedef void tFnProcessf_blk(void *pScr, void *y, const void *R, const void *D, const void *A, const void *B);
    int32_t M, N, P, L, Nq;
    int qY, qB, qR, qA, qYBRA;
    void *A, *R, *D, *B, *Y;
    void *w;// reciprocal of main diagonal in floating point representation
    const tCholApi *pCholContext;
    te_fun_ptr_t fxn;
    fxn = context->target.fut;

    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    P = MAX(0,context->args.dim[2]);
    L = MAX(0,context->args.dim[3]);
    (void)M, (void)N, (void)P, (void)L;
    R = vecGetElem(&context->dataSet.X, 0);
    A = vecGetElem(&context->dataSet.U, 0);
    B = vecGetElem(&context->dataSet.V, 0);
    Y = vecGetElem(&context->dataSet.Z, 0);
    D = vecGetElem(&context->dataSet.auxVec[1], 0);
    w = vecGetElem(&context->dataSet.Y,0);

    qYBRA = 0;
    Nq = te_szHelper_cholfwd(context).NQ;
    if (Nq > 0)
    {
        int32_t *Q;
        Q = (( int32_t *)vecGetElem(&context->dataSet.auxVec[0], 0));
        if (Nq == 4)/* qA, qB, qX, qY */
        {
            qA=Q[0];
            qB=Q[1];
            qY=Q[3];
            qYBRA = qY - qB;
        }
        else/* qA, qB, qR, qY, qX */
        {
            qA=Q[0];
            qB=Q[1];
            qR=Q[2];
            qY=Q[3];
            qYBRA = qY-qB+qR-qA;
        }
    }
    pCholContext = (const tCholApi *)context->desc->extraPtr;
    /* Copy w to D with right formatting */
    pCholContext->te_cholFxn_unpackD(D, w, N, L);

    /* Apply function */
    {
        tVec vScr;
        int scr_sz;

        /* Allocate scratch buffer */
        scr_sz = pCholContext->fwdScratchSz();
        if (!vecAlloc(&vScr, scr_sz, context->desc->isAligned, FMT_REAL|FMT_UINT8, NULL)) 
        {
            printf( "te_processFxn_cholfwd(): failed to allocate memory\n");
        }
        else
        {
            void * pScr;
            pScr = vecGetElem(&vScr, 0);
            {    /* reporting */
                tReportFUT fut[2];
                fut[0]=(tReportFUT)context->target.fut;
                fut[1]=(tReportFUT)pCholContext->fwdScratchSz;
                vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
            }
            switch(context->desc->extraParam)
            {
                case CHOLAPI_FIX_SZ|CHOLAPI_PACKED:
                    if ((context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT32 || (context->desc->fmt & FMT_DTYPE_MASK) == FMT_FLOAT16)
                        ((tFnProcessf_blk*)fxn)(pScr,Y,R,D,A,B);
                    else
                        ((tFnProcess*)fxn)(pScr,Y,R,D,A,B,qYBRA);                        
                    break;
                default: ASSERT(0);
            }
            vecFree(&vScr);
        }
    }
} /* te_processFxn_cholfwd() */

/* Allocate vectors and load the data set for Cholesky MMSE functions */
int te_loadFxn_cholmmse( tTestEngContext * context )
{
    int M,N,P,L,Narg,res = 0,Nq=0,NTotalArgs;
    int Sa,Sb,Sx;
    int sz=0;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    P = MAX(0,context->args.dim[2]);
    L = MAX(0,context->args.dim[3]);

    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
        case FMT_FLOAT32: sz=sizeof(float32_t); Nq=0; break;
        case FMT_FRACT32: sz=sizeof(fract32);   Nq=5; /* qA, qB, qR, qY, qX*/
        break;
    default: ASSERT(0);
    }

    /* Allocate input data vector memory. */
    if ((context->desc->extraParam & CHOLAPI_LAYOUT_MASK) == CHOLAPI_PACKED)
    {
        Sx=(context->desc->fmt & FMT_CPLX) ? getSpace(N*P, 2*sz) : getSpace(N*P, sz);
        Sa=(context->desc->fmt & FMT_CPLX) ? getSpace(M*N, 2*sz) : getSpace(M*N, sz);
        Sb=(context->desc->fmt & FMT_CPLX) ? getSpace(M*P, 2*sz) : getSpace(M*P, sz);
    }
    else
    {
        Sx=N*P;
        Sa=M*N;
        Sb=M*P;
    }

    Narg = 0;
    Narg+=vecsAlloc( context->desc->isAligned, (context->desc->fmt & FMT_DTYPE_MASK),
                    &context->dataSet.U  , L,
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X  , Sa*L,    // A
                    &context->dataSet.Y  , Sb*L,    // B
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.Z  , Sx*L,    // X
                    &context->dataSet.Zlo, Sx*L,    // X
                    &context->dataSet.Zhi, Sx*L,    // X
                    0);
    NTotalArgs=6;
    /* allocate memory for parameters Q */
    {
        Narg+=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.auxVec[0]  , Nq, /* Q */
                        0);
        NTotalArgs += 1;
    }
    if(Narg!=NTotalArgs)
    {
        printf( "te_loadFxn_cholmmses(): failed to allocate memory\n" );
    }
    else
    {
        /* clean output of function just to protect from ferret warnings appearing 
           when test environment reads data in the holes between matrices */
        memset(vecGetElem(&context->dataSet.Z,0),0,vecGetSize(&context->dataSet.Z));

        /* Load vectors data from the SEQ-file. */
        if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.U,
                             &context->dataSet.auxVec[0],
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 
                             0 ) )
        {
            printf( "te_loadFxn_cholmmses(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholmmse() */

/* processing function for Cholesky MMSE functions */
void te_processFxn_cholmmse(tTestEngContext * context)
{
    typedef void tFnProcess(
                void * pScr,
                void * restrict x,
          const void * restrict A,
          const void * restrict B,
                  float32_t sigma2);
    typedef void tFnProcess32x32(
                void * pScr,
                void * restrict x,
          const void * restrict A,
          const void * restrict B,
                int32_t sigma2, int qRA,int qYBRA,int qXYR);
    void* pScr;
    const tCholApi * pApi=(const tCholApi*)context->desc->extraPtr;
    size_t scrsz;

    int32_t N,M,P,L;
    void *A,*B,*X;
    void *sigma2;
    tFnProcess *fxn;
    const int32_t *q;
    fxn=(tFnProcess *)context->target.fut;

    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    P = MAX(0,context->args.dim[2]);
    L = MAX(0,context->args.dim[3]);
    (void)M,(void)N,(void)P,(void)L;
    scrsz = pApi->mmseScratchSz();
    pScr = mallocAlign(scrsz,0);
    if (pScr==0) 
    {
        printf( "te_processFxn_cholmmsesf(): failed to allocate memory\n");
        return;
    }

    A      =vecGetElem(&context->dataSet.X  ,0);
    B      =vecGetElem(&context->dataSet.Y  ,0);
    sigma2 =vecGetElem(&context->dataSet.U  ,0);
    X      =vecGetElem(&context->dataSet.Z  ,0);
    q      =(const int32_t*)vecGetElem(&context->dataSet.auxVec[0],0);


    {    /* reporting */
        tReportFUT fut[2];
        fut[0]=(tReportFUT)context->target.fut;
        fut[1]=(tReportFUT)pApi->mmseScratchSz;
        vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
    }
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    { 
        case FMT_FLOAT32:
            fxn(pScr,X,A,B,((const float32_t*)sigma2)[0]); break;
        case FMT_FRACT32:
            ((tFnProcess32x32*)fxn)(pScr,X,A,B,
                 ((const int32_t*)sigma2)[0],
                 q[2]-q[0],
                 q[3]-q[1]+q[2]-q[0],
                 q[4]-q[3]+q[2]); break;
        default: ASSERT(0);
    }
    freeAlign(pScr);
} /* te_processFxn_cholmmse() */

/* Allocate vectors and load the data set for cholnf functions */
int te_loadFxn_cholf( tTestEngContext * context )
{
    int M,N,Nr,Nd,Narg,res = 0;
    int Sa;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);

    /* Allocate input data vector memory. */
    Nr=(context->desc->fmt & FMT_CPLX) ? getSpace(N*(N+1),sizeof(complex_float))>>1   : getSpace(N*(N+1)>>1,sizeof(float32_t));
    Nd=(context->desc->fmt & FMT_CPLX) ? getSpace(N  *2  ,sizeof(complex_float))>>1   : getSpace(N         ,sizeof(float32_t));
    Sa=(context->desc->fmt & FMT_CPLX) ? getSpace(M*N*2  ,sizeof(complex_float))>>1   : getSpace(M*N       ,sizeof(float32_t));
    Narg =vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X  , Sa,
                    &context->dataSet.Z  , Nr,
                    &context->dataSet.Zlo, Nr,
                    &context->dataSet.Zhi, Nr,
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, FMT_FLOAT32,
                    &context->dataSet.Y  , 1,   /* sigma2 */
                    &context->dataSet.W  , Nd,
                    &context->dataSet.Wlo, Nd,
                    &context->dataSet.Whi, Nd,
                    0 );

    if(Narg!=8)
    {
        printf( "te_loadFxn_cholf(): failed to allocate memory\n" );
    }
    else
    {
        /* Load vectors data from the SEQ-file. */
        if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.Y,
                             &context->dataSet.X,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 
                             &context->dataSet.Wlo,
                             &context->dataSet.Whi, 
                             0 ) )
        {
            printf( "te_loadFxn_cholf(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholf() */


/* Allocate vectors and load the data set for cholfwdf functions */
int te_loadFxn_cholfwdf( tTestEngContext * context )
{
    int M,N,P,Nr,Nd,Narg,res = 0;
    int Sa,Sb,Sy;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    P = MAX(0,context->args.dim[2]);

    /* Allocate input data vector memory. */
    Nr=(context->desc->fmt & FMT_CPLX) ? getSpace((N*(N+1)),sizeof(complex_float))>>1 : getSpace((N*(N+1)>>1),sizeof(float32_t));
    Nd=(context->desc->fmt & FMT_CPLX) ? getSpace(N  *2    ,sizeof(complex_float))>>1 : getSpace(N           ,sizeof(float32_t));
    Sy=(context->desc->fmt & FMT_CPLX) ? getSpace(N*P*2    ,sizeof(complex_float))>>1 : getSpace(N*P         ,sizeof(float32_t));
    Sa=(context->desc->fmt & FMT_CPLX) ? getSpace(M*N*2    ,sizeof(complex_float))>>1 : getSpace(M*N         ,sizeof(float32_t));
    Sb=(context->desc->fmt & FMT_CPLX) ? getSpace(M*P*2    ,sizeof(complex_float))>>1 : getSpace(M*P         ,sizeof(float32_t));

    Narg =0;
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X  , Nr,    // R
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, FMT_FLOAT32,
                    &context->dataSet.Y  , Nd,
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.U  , Sa,    // A
                    &context->dataSet.V  , Sb,    // B
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.Z   , Sy,    // Y
                    &context->dataSet.Zlo , Sy,    // Y
                    &context->dataSet.Zhi , Sy,    // Y
                    0);

    if(Narg!=7)
    {
        printf( "te_loadFxn_cholfwdf(): failed to allocate memory\n" );
    }
    else
    {
        /* Load vectors data from the SEQ-file. */
        if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.U,
                             &context->dataSet.V,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 
                             0 ) )
        {
            printf( "te_loadFxn_cholfwdf(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholfwdf() */


/* Allocate vectors and load the data set for cholbkwf functions */
int te_loadFxn_cholbkwf( tTestEngContext * context )
{
    int N,P,Nr,Nd,Narg,res = 0;
    int Sx,Sy;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    N = MAX(0,context->args.dim[0]);
    P = MAX(0,context->args.dim[1]);

    /* Allocate input data vector memory. */
    Nr=(context->desc->fmt & FMT_CPLX) ? getSpace((N*(N+1)),sizeof(complex_float))>>1 :getSpace((N*(N+1))>>1,sizeof(float32_t));
    Nd=(context->desc->fmt & FMT_CPLX) ? getSpace(N  *2    ,sizeof(complex_float))>>1 :getSpace(N           ,sizeof(float32_t));
    Sy=(context->desc->fmt & FMT_CPLX) ? getSpace(N*P*2    ,sizeof(complex_float))>>1 :getSpace(N*P         ,sizeof(float32_t));
    Sx=(context->desc->fmt & FMT_CPLX) ? getSpace(N*P*2    ,sizeof(complex_float))>>1 :getSpace(N*P         ,sizeof(float32_t));
    Narg =0;
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.X  , Nr,    // R
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, FMT_FLOAT32,
                    &context->dataSet.Y  , Nd,
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.U  , Sy,    // Y
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.Z   , Sx,    // X
                    &context->dataSet.Zlo , Sx,    // X
                    &context->dataSet.Zhi , Sx,    // X
                    0);

    if(Narg!=6)
    {
        printf( "te_loadFxn_cholbkwf(): failed to allocate memory\n" );
    }
    else
    {
        /* Load vectors data from the SEQ-file. */
        if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.U,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 
                             0 ) )
        {
            printf( "te_loadFxn_cholbkwf(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholbkwf() */


/* Allocate vectors and load the data set for cholmmsef functions */
int te_loadFxn_cholmmsef( tTestEngContext * context )
{
    int M,N,P,Narg,res = 0;
    int Sa,Sb,Sx;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    M = MAX(0,context->args.dim[0]);
    N = MAX(0,context->args.dim[1]);
    P = MAX(0,context->args.dim[2]);

    /* Allocate input data vector memory. */
    Sx=(context->desc->fmt & FMT_CPLX) ? getSpace(N*P*2,sizeof(complex_float))>>1 : getSpace(N*P,sizeof(float32_t));
    Sa=(context->desc->fmt & FMT_CPLX) ? getSpace(M*N*2,sizeof(complex_float))>>1 : getSpace(M*N,sizeof(float32_t));
    Sb=(context->desc->fmt & FMT_CPLX) ? getSpace(M*P*2,sizeof(complex_float))>>1 : getSpace(M*P,sizeof(float32_t));

    Narg =0;
    Narg+=vecsAlloc( context->desc->isAligned, FMT_FLOAT32,
                    &context->dataSet.Y  , 1,    // sigma2
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.U  , Sa,    // A
                    &context->dataSet.V  , Sb,    // B
                    0);
    Narg+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                    &context->dataSet.Z   , Sx,    // x
                    &context->dataSet.Zlo , Sx,    // x
                    &context->dataSet.Zhi , Sx,    // x
                    0);

    if(Narg!=6)
    {
        printf( "te_loadFxn_cholmmsef(): failed to allocate memory\n" );
    }
    else
    {
        /* Load vectors data from the SEQ-file. */
        if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.U,
                             &context->dataSet.V,
                             &context->dataSet.Y,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 
                             0 ) )
        {
            printf( "te_loadFxn_cholmmsef(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )    te_freeVectors(context);  
    return (res);
} /* te_loadFxn_cholmmsef() */

/* Allocate vectors and load the data set for cholpreprocessf functions */
int te_loadFxn_cholpreprocessf(tTestEngContext * context)
{
    int M, N, Nr=0, Narg, res = 0, Nq=0,fmtA=0,fmtZ=0,fmtSigma2=0;
    int Sa=0;
    ASSERT(context && context->seqFile);
    memset(&context->dataSet, 0, sizeof(context->dataSet));
    M = MAX(0, context->args.dim[0]);
    N = MAX(0, context->args.dim[1]);

    /* Allocate input data vector memory. */
    switch(context->desc->fmt & (FMT_DTYPE_MASK|FMT_CPLX|FMT_REAL))
    {
    case FMT_FLOAT32|FMT_REAL:
        Nr = getSpace(N*N,sizeof(float32_t));
        Sa = getSpace(M*N,sizeof(float32_t));
        Nq = 0;
        fmtA=fmtZ= FMT_FLOAT32|FMT_REAL;
        fmtSigma2=FMT_FLOAT32|FMT_REAL;
        break;
    case FMT_FLOAT32|FMT_CPLX:
        Nr =  getSpace(N*N * 2,sizeof(complex_float)) >> 1;
        Sa =  getSpace(M*N * 2,sizeof(complex_float)) >> 1;
        Nq = 0;
        fmtA=fmtZ= FMT_FLOAT32|FMT_CPLX;
        fmtSigma2=FMT_FLOAT32|FMT_REAL;
        break;
    case FMT_FRACT32|FMT_REAL:
        Nr = getSpace(N*N,sizeof(int32_t));
        Sa = getSpace(M*N,sizeof(int32_t));
        Nq = 5;
        fmtA =FMT_FRACT32|FMT_REAL;
        fmtZ =FMT_INT64  |FMT_REAL;
        fmtSigma2=FMT_FRACT32|FMT_REAL;
        break;
    case FMT_FRACT32|FMT_CPLX:
        Nr =  getSpace(N*N * 2,sizeof(complex_fract64)) >> 1;
        Sa =  getSpace(M*N * 2,sizeof(complex_fract32)) >> 1;
        Nq = 5;
        fmtA =FMT_FRACT32|FMT_CPLX;
        fmtZ =FMT_INT64  |FMT_CPLX;
        fmtSigma2=FMT_FRACT32|FMT_REAL;
        break;
        default: ASSERT(0);
    }
    Narg=0;
    Narg+= vecsAlloc(context->desc->isAligned, fmtA, &context->dataSet.X, Sa, NULL);
    Narg+= vecsAlloc(context->desc->isAligned, FMT_INT32, &context->dataSet.auxVec[0], Nq, NULL);
    Narg+= vecsAlloc(context->desc->isAligned, fmtZ, 
        &context->dataSet.Z, Nr,
        &context->dataSet.Zlo, Nr,
        &context->dataSet.Zhi, Nr,
        NULL);
    Narg += vecsAlloc(context->desc->isAligned, fmtSigma2, &context->dataSet.Y, 1,NULL);   /* sigma2 */
    if (Narg != 6)
    {
        printf("te_loadFxn_cholpreprocessf(): failed to allocate memory\n");
    }
    else
    {
        /* Load vectors data from the SEQ-file. */
        if (!seqFileReadVecs(context->seqFile,
            &context->dataSet.X,
            &context->dataSet.Y,
            &context->dataSet.auxVec[0],
            &context->dataSet.Zlo,
            &context->dataSet.Zhi,
            0))
        {
            printf("te_loadFxn_cholpreprocessf(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if (!res)    te_freeVectors(context);
    return (res);
} /* te_loadFxn_cholpreprocessf() */

/* Allocate vectors and load the data set for cholpinvf functions */
int te_loadFxn_cholpinvf(tTestEngContext * context)
{
    int M, N, Nq=0, fmtA=0, fmtX=0, fmtSigma2=0, Narg, res = 0;
    int Sa=0,Nx=0;
    ASSERT(context && context->seqFile);
    memset(&context->dataSet, 0, sizeof(context->dataSet));
    M = MAX(0, context->args.dim[0]);
    N = MAX(0, context->args.dim[1]);
    switch(context->desc->fmt & (FMT_DTYPE_MASK|FMT_CPLX|FMT_REAL))
    {
    case FMT_FLOAT32|FMT_REAL:
        Nx = getSpace(N*N,sizeof(float32_t));
        Sa = getSpace(M*N,sizeof(float32_t));
        Nq = 0;
        fmtA=fmtX= FMT_FLOAT32|FMT_REAL;
        fmtSigma2=FMT_FLOAT32|FMT_REAL;
        break;
    case FMT_FLOAT32|FMT_CPLX:
        Nx =  getSpace(N*N * 2,sizeof(complex_float)) >> 1;
        Sa =  getSpace(M*N * 2,sizeof(complex_float)) >> 1;
        Nq = 0;
        fmtA=fmtX= FMT_FLOAT32|FMT_CPLX;
        fmtSigma2=FMT_FLOAT32|FMT_REAL;
        break;
    case FMT_FRACT32|FMT_REAL:
        Nx = getSpace(N*N,sizeof(int32_t));
        Sa = getSpace(M*N,sizeof(int32_t));
        Nq = 5;
        fmtA =FMT_FRACT32|FMT_REAL;
        fmtX =FMT_FRACT32|FMT_REAL;
        fmtSigma2=FMT_FRACT32|FMT_REAL;
        break;
    case FMT_FRACT32|FMT_CPLX:
        Nx =  getSpace(N*N * 2,sizeof(complex_fract64)) >> 1;
        Sa =  getSpace(M*N * 2,sizeof(complex_fract32)) >> 1;
        Nq = 5;
        fmtA =FMT_FRACT32|FMT_CPLX;
        fmtX =FMT_FRACT32|FMT_REAL;
        fmtSigma2=FMT_FRACT32|FMT_REAL;
        break;
        default: ASSERT(0);
    }
    /* Allocate input data vector memory. */
    Narg=0;
    Narg+= vecsAlloc(context->desc->isAligned, fmtA,&context->dataSet.X, Sa,NULL);
    Narg+= vecsAlloc(context->desc->isAligned, fmtSigma2,&context->dataSet.Y, 1,NULL);   /* sigma2 */
    Narg+= vecsAlloc(context->desc->isAligned, FMT_INT32,&context->dataSet.auxVec[0], Nq,NULL); // Q
    Narg+= vecsAlloc(context->desc->isAligned, context->desc->fmt,
        &context->dataSet.Z, Nx,
        &context->dataSet.Zlo, Nx,
        &context->dataSet.Zhi, Nx,
        NULL);
    if (Narg != 6)
    {
        printf("te_loadFxn_cholpinvf(): failed to allocate memory\n");
    }
    else
    {
        /* Load vectors data from the SEQ-file. */
        if (!seqFileReadVecs(context->seqFile,
        &context->dataSet.X,
        &context->dataSet.Y,
        &context->dataSet.auxVec[0],
        &context->dataSet.Zlo,
        &context->dataSet.Zhi,
        0))
        {
            printf("te_loadFxn_cholpinvf(): failed to read vectors data\n");
        }
        else
        {
            res = 1;
        }
    }
    /* Free vectors data if failed. */
    if (!res)    te_freeVectors(context);
    return (res);
} /* te_loadFxn_cholpinvf() */

/* processing function for cholf */
void te_processFxn_cholf( tTestEngContext * context )
{
    typedef void tFnProcessmxn(
            void * pScr,
            complex_float * restrict R, 
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t sigma2);
    //tAllocPtr scratch;
    tVec vD;
    void *pScr;
    int32_t N,Nd;
    int nVec;
    complex_float *A,*R,*D;
    float32_t sigma2;
    const tCholApi* pApi=(const tCholApi*)context->desc->extraPtr;
    float32_t *w;   // reciprocal of main diagonal in floating point representation
    tFnProcessmxn *fxn = (tFnProcessmxn *)context->target.fut;

    N = MAX(0,context->args.dim[1]);
    Nd = (context->desc->fmt & FMT_CPLX) ? getSpace(N*2,sizeof(complex_float))>>1 : getSpace(N,sizeof(complex_float));
    pScr = mallocAlign(pApi->cholScratchSz(),1<<(LOG2_SIMD_WIDTH+1));
    nVec = vecAlloc( &vD, Nd, context->desc->isAligned, context->desc->fmt, NULL);
    if (pScr==NULL || nVec==0)  
    {
        printf( "te_processFxn_cholf(): failed to allocate memory\n");
        return;
    }

    A     =((complex_float*)vecGetElem(&context->dataSet.X  ,0));
    sigma2=(( float32_t   *)vecGetElem(&context->dataSet.Y  ,0))[0];
    R     =((complex_float*)vecGetElem(&context->dataSet.Z  ,0));
    w     =(( float32_t *)vecGetElem(&context->dataSet.W  ,0));
    D     =((complex_float*)vecGetElem(&vD  ,0));
    memset(D,0,vecGetSize(&vD));
	{	/* reporting */
		tReportFUT fut[2];
		fut[0]=(tReportFUT)context->target.fut;
		fut[1]=(tReportFUT)pApi->cholScratchSz;
		vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
	}

    fxn(pScr, R, D, A, sigma2);

    /* copy D to w with right formatting */
    if (context->desc->fmt & FMT_CPLX)
    {
        int n;
        for (n = 0; n < Nd; n++) w[n] = ((float32_t*)D)[2 * n];
    }
    else
    {
        int n;
        for (n = 0; n < Nd; n++) w[n] = ((float32_t*)D)[n];
    }

    freeAlign(pScr);
    vecFree(&vD);
} /* te_processFxn_cholf() */


/* processing function for cholfwdf */
void te_processFxn_cholfwdf( tTestEngContext * context )
{
    typedef void tFnProcessmxn(
            void *pScr,
                  float32_t* restrict y,
            const float32_t* restrict R, 
            const float32_t* restrict D,
            const float32_t* restrict A, 
            const float32_t* restrict B);
    const tCholApi* pApi=(const tCholApi*)context->desc->extraPtr;
    //tAllocPtr scratch;
    tVec vD;
    void *pScr;
    int32_t N,Nd;
    int nVec;
    float32_t *A,*R,*D,*B,*Y;
    float32_t *w;   // reciprocal of main diagonal in floating point representation
    tFnProcessmxn *fxn = (tFnProcessmxn *)context->target.fut;

    N = MAX(0,context->args.dim[1]);

    pScr= mallocAlign(/*&scratch,*/pApi->fwdScratchSz(),1<<(LOG2_SIMD_WIDTH+1));
    Nd = (context->desc->fmt & FMT_CPLX) ? getSpace(N*2,sizeof(complex_float))>>1 : getSpace(N,sizeof(float32_t));
    nVec=vecsAlloc( context->desc->isAligned, context->desc->fmt, &vD, Nd, NULL);
    if (pScr==NULL || nVec==0) 
    {
        printf( "te_processFxn_cholfwdf(): failed to allocate memory\n");
        return;
    }

    R =((float32_t*)vecGetElem(&context->dataSet.X  ,0));
    w =((float32_t*)vecGetElem(&context->dataSet.Y  ,0));
    A =((float32_t*)vecGetElem(&context->dataSet.U  ,0));
    B =((float32_t*)vecGetElem(&context->dataSet.V  ,0));
    Y =((float32_t*)vecGetElem(&context->dataSet.Z  ,0));
    D =((float32_t*)vecGetElem(&vD  ,0));
    /* copy w to D with right formatting */
    if (context->desc->fmt & FMT_CPLX)
    {
        int n;
        for (n = 0; n < N; n++) D[2 * n + 1] = D[2 * n + 0] = w[n];
    }
    else
    {
        int n;
        for (n = 0; n < N; n++) D[n] = w[n];
    }
	{	/* reporting */
		tReportFUT fut[2];
		fut[0]=(tReportFUT)context->target.fut;
		fut[1]=(tReportFUT)pApi->fwdScratchSz;
		vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
	}

    fxn(pScr, Y, R, D, A, B);

    freeAlign(pScr);
    vecFree(&vD);
} /* te_processFxn_cholfwdf() */


/* processing function for cholbkwf */
void te_processFxn_cholbkwf( tTestEngContext * context )
{
    typedef void tFnProcessmxn(
            void *pScr,
                  float32_t* restrict x, 
            const float32_t* restrict R,
            const float32_t* restrict D,
            const float32_t* restrict y);
    const tCholApi* pApi=(const tCholApi*)context->desc->extraPtr;
    //tAllocPtr scratch;
    tVec vD;
    void *pScr;
    int32_t N,Nd;
    int nVec;
    float32_t *R,*D,*X,*Y;
    float32_t *w;   // reciprocal of main diagonal in floating point representation
    tFnProcessmxn *fxn = (tFnProcessmxn *)context->target.fut;

    N = MAX(0,context->args.dim[0]);

    pScr= mallocAlign(/*&scratch,*/pApi->bkwScratchSz(),1<<(LOG2_SIMD_WIDTH+1));
    Nd = (context->desc->fmt & FMT_CPLX) ? getSpace(N*2,sizeof(complex_float))>>1 : getSpace(N,sizeof(float32_t));
    nVec=vecsAlloc( context->desc->isAligned, context->desc->fmt, &vD, Nd, NULL);
    if (pScr==NULL || nVec==0) 
    {
        printf( "te_processFxn_cholbkwf(): failed to allocate memory\n");
        return;
    }

    R =((float32_t*)vecGetElem(&context->dataSet.X  ,0));
    w =((float32_t*)vecGetElem(&context->dataSet.Y  ,0));
    Y =((float32_t*)vecGetElem(&context->dataSet.U  ,0));
    X =((float32_t*)vecGetElem(&context->dataSet.Z  ,0));
    D =((float32_t*)vecGetElem(&vD  ,0));
    /* copy w to D with right formatting */
    if (context->desc->fmt & FMT_CPLX)
    {
        int n;
        for (n = 0; n < N; n++) D[2 * n + 1] = D[2 * n + 0] = w[n];
    }
    else
    {
        int n;
        for (n = 0; n < N; n++) D[n] = w[n];
    }
	{	/* reporting */
		tReportFUT fut[2];
		fut[0]=(tReportFUT)context->target.fut;
		fut[1]=(tReportFUT)pApi->bkwScratchSz;
		vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
	}

    fxn(pScr, X, R, D, Y);

    freeAlign(pScr);
    vecFree(&vD);
} /* te_processFxn_cholbkwf() */


/* processing function for cholmmsef */
void te_processFxn_cholmmsef( tTestEngContext * context )
{
    typedef void tFnProcessmxn(
            void *pScr,
                  float32_t* restrict x,
            const float32_t* restrict A, 
            const float32_t* restrict B, 
            const float32_t sigma2);
    const tCholApi* pApi=(const tCholApi*)context->desc->extraPtr;
    //tAllocPtr scratch;
    void *pScr;
    float32_t *A,*B,*X;
    float32_t sigma2;
    tFnProcessmxn *fxn = (tFnProcessmxn *)context->target.fut;

    pScr= mallocAlign(/*&scratch,*/pApi->mmseScratchSz(),1<<(LOG2_SIMD_WIDTH+1));
    if (pScr==NULL ) 
    {
        printf( "te_processFxn_cholmmsef(): failed to allocate memory\n");
        return;
    }
    sigma2 =((float32_t*)vecGetElem(&context->dataSet.Y  ,0))[0];
    A      =((float32_t*)vecGetElem(&context->dataSet.U  ,0));
    B      =((float32_t*)vecGetElem(&context->dataSet.V  ,0));
    X      =((float32_t*)vecGetElem(&context->dataSet.Z  ,0));
    {	/* reporting */
        tReportFUT fut[2];
        fut[0]=(tReportFUT)context->target.fut;
        fut[1]=(tReportFUT)pApi->mmseScratchSz;
        vReportAdd(fut,2,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
    }

    fxn(pScr, X, A, B, sigma2);

    freeAlign(pScr);
} /* te_processFxn_cholmmsef() */

/* processing function for cholpreprocessf */
void te_processFxn_cholpreprocessf(tTestEngContext * context)
{
    typedef void tFnProcessmxn(void * pScr,void * restrict R,const void * restrict A,float32_t sigma2);
    typedef void tFnProcessmxn_fixpoint(void * pScr,void * restrict R,const void* restrict A,int32_t sigma2, int qRA);
    void *pScr;
    void *A, *R;
    const void* sigma2;
    const int32_t *q;
    const tCholApi* pApi = (const tCholApi*)context->desc->extraPtr;
    tFnProcessmxn *fxn = (tFnProcessmxn *)context->target.fut;

    pScr = mallocAlign(pApi->preprocessScratchSz(), 1 << (LOG2_SIMD_WIDTH + 1));
    if (pScr == NULL)
    {
        printf("te_processFxn_cholpreprocessf(): failed to allocate memory\n");
        return;
    }

    A = vecGetElem(&context->dataSet.X, 0);
    sigma2 = vecGetElem(&context->dataSet.Y, 0);
    R = vecGetElem(&context->dataSet.Z, 0);
    q = (const int32_t*)vecGetElem(&context->dataSet.auxVec[0], 0);
    {	/* reporting */
        tReportFUT fut[2];
        fut[0] = (tReportFUT)context->target.fut;
        fut[1] = (tReportFUT)pApi->preprocessScratchSz;
        vReportAdd(fut, 2, NULL, context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
    }

    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
        case FMT_FLOAT32:
            fxn(pScr, R, A, ((const float32_t*)sigma2)[0]);
            break;
        case FMT_FRACT32:
            ((tFnProcessmxn_fixpoint*)fxn)(pScr, R, A, ((const int32_t*)sigma2)[0],q[2]-q[0]);
            break;
        default: ASSERT(0);
    }
    freeAlign(pScr);
} /* te_processFxn_cholpreprocessf() */

/* processing function for cholpinvf */
void te_processFxn_cholpinvf(tTestEngContext * context)
{
    typedef void tFnProcessmxn(
        void * pScr,
        void * restrict x,
        const void* restrict A,
        float32_t sigma2);
    typedef void tFnProcessmxn_fixedpoint(
        void * pScr,
        void * restrict x,
        const void * restrict A,
        const int32_t sigma2,int qRA, int qYBRA, int qXYR);
    void *pScr;
    void *A, *x;
    void *sigma2;
    int32_t *q;
    const tCholApi* pApi = (const tCholApi*)context->desc->extraPtr;
    tFnProcessmxn *fxn = (tFnProcessmxn *)context->target.fut;

    pScr = mallocAlign(pApi->pinvScratchSz(), 1 << (LOG2_SIMD_WIDTH + 1));
    if (pScr == NULL)
    {
        printf("te_processFxn_cholpinvf(): failed to allocate memory\n");
        return;
    }

    A = vecGetElem(&context->dataSet.X, 0);
    x = vecGetElem(&context->dataSet.Z, 0);
    sigma2 = vecGetElem(&context->dataSet.Y, 0);
    q      = (int32_t*)vecGetElem(&context->dataSet.auxVec[0], 0);
    {	/* reporting */
        tReportFUT fut[2];
        fut[0] = (tReportFUT)context->target.fut;
        fut[1] = (tReportFUT)pApi->pinvScratchSz;
        vReportAdd(fut, 2, NULL, context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
    }
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
        case FMT_FLOAT32:
            fxn(pScr, x, A, ((float32_t*)sigma2)[0]);
            break;
        case FMT_FRACT32:
            ((tFnProcessmxn_fixedpoint*)fxn)(pScr, x, A, ((int32_t*)sigma2)[0],
                q[2]-q[0], q[3]-q[1] + q[2]-q[0], q[4]-q[3]+q[2]);
            break;
    }
    freeAlign(pScr);
} /* te_processFxn_cholpinvf() */


#define MAX_FUNC_NUM 4
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extraParam, extraPtr,argNum, align, loadFxn, procFxn ) { (fmt),extraParam,extraPtr,(argNum),(align),NULL,NULL,(loadFxn),(procFxn) }

static const struct 
{
  tTestEngTarget funcList[MAX_FUNC_NUM];
  tTestEngDesc   testDesc;
}
testDefTbl[] =
{
    // 32-bit fixed point
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp4x4_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp6x6_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp8x8_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp10x10_32x32     ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp4x4_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp6x6_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp8x8_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp10x10_32x32      ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },

    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst4x4_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol4x4_32x32_Api  , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst6x6_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol6x6_32x32_Api  , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst8x8_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol8x8_32x32_Api  , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst10x10_32x32 ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst4x4_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol4x4_32x32_Api   , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst6x6_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol6x6_32x32_Api   , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst8x8_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol8x8_32x32_Api   , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst10x10_32x32  ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol10x10_32x32_Api , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },

    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst4x4_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol4x4_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst6x6_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol6x6_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst8x8_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol8x8_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst10x10_32x32 ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol10x10_32x32_Api, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst4x4_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol4x4_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst6x6_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol6x6_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst8x8_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol8x8_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst10x10_32x32  ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol10x10_32x32_Api , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },

    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver4x4_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver6x6_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver8x8_32x32   ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver10x10_32x32 ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver4x4_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver6x6_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver8x8_32x32    ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver10x10_32x32  ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },

    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess4x4_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess6x6_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess8x8_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess10x10_32x32     ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess4x4_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess6x6_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess8x8_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess10x10_32x32      ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api,  TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },

    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv4x4_32x32    ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv6x6_32x32    ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv8x8_32x32    ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv10x10_32x32  ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv4x4_32x32     ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv6x6_32x32     ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv8x8_32x32     ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv10x10_32x32   ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api,  TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },

    // floating point
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp4x4f            ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp6x6f            ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp8x8f            ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp10x10f          ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst4x4f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst6x6f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst8x8f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst10x10f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst4x4f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst6x6f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst8x8f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst10x10f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver4x4f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver6x6f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver8x8f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver10x10f      ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess4x4f  ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess6x6f  ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess8x8f  ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess10x10f), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv4x4f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv6x6f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv8x8f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv10x10f      ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
	
    { FUNC_LIST((tTestEngTarget)&matcholdecomp4x4f            ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp6x6f            ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp8x8f            ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp10x10f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst4x4f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst6x6f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst8x8f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst10x10f       ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst4x4f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst6x6f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst8x8f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst10x10f       ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver4x4f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver6x6f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver8x8f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver10x10f      ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess4x4f  ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess6x6f  ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess8x8f  ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess10x10f), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv4x4f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv6x6f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv8x8f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv10x10f      ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },

    { FUNC_LIST(NULL), TEST_DESC(0,0,0,0,0, NULL, NULL) } /* End of table */
};

int main_chol  ( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int res = 1;
#define TEST_SINGLE 0
#define TEST_DOUBLE 1
#define DO_TEST(fxn, seqFile)                                                                    \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), (seqFile),         \
                                                       isFull, isVerbose, breakOnError ) )

  if (phaseNum==0 || phaseNum==1)
  {
        DO_TEST( &cmatcholdecomp4x4_32x32          , "cmatcholdecomp4x4_32x32_cond20dynRange29qRAm3_0.seq"      );
        DO_TEST( &cmatcholdecomp4x4_32x32          , "cmatcholdecomp4x4_32x32_cond200dynRange29qRAm3_0.seq"     );
        DO_TEST( &cmatcholdecomp4x4_32x32          , "cmatcholdecomp4x4_32x32_cond2000dynRange29qRAm3_0.seq"    );
        DO_TEST( &cmatcholdecomp6x6_32x32          , "cmatcholdecomp6x6_32x32_cond20dynRange29qRAm3_0.seq"      );
        DO_TEST( &cmatcholdecomp6x6_32x32          , "cmatcholdecomp6x6_32x32_cond200dynRange29qRAm3_0.seq"     );
        DO_TEST( &cmatcholdecomp6x6_32x32          , "cmatcholdecomp6x6_32x32_cond2000dynRange29qRAm3_0.seq"    );
        DO_TEST( &cmatcholdecomp8x8_32x32          , "cmatcholdecomp8x8_32x32_cond20dynRange29qRAm3_0.seq"      );
        DO_TEST( &cmatcholdecomp8x8_32x32          , "cmatcholdecomp8x8_32x32_cond200dynRange29qRAm3_0.seq"     );
        DO_TEST( &cmatcholdecomp8x8_32x32          , "cmatcholdecomp8x8_32x32_cond2000dynRange29qRAm3_0.seq"    );
        DO_TEST( &cmatcholdecomp10x10_32x32        , "cmatcholdecomp10x10_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &cmatcholdecomp10x10_32x32        , "cmatcholdecomp10x10_32x32_cond200dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholdecomp10x10_32x32        , "cmatcholdecomp10x10_32x32_cond2000dynRange29qRAm3_0.seq"  );
        DO_TEST( &cmatcholfwdsubst4x4_32x32        , "cmatcholfwdsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholfwdsubst6x6_32x32        , "cmatcholfwdsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholfwdsubst8x8_32x32        , "cmatcholfwdsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholfwdsubst10x10_32x32      , "cmatcholfwdsubst10x10_32x32_cond20dynRange29qRAm3_0.seq" );
        DO_TEST( &cmatcholbkwsubst4x4_32x32        , "cmatcholbkwsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholbkwsubst6x6_32x32        , "cmatcholbkwsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholbkwsubst8x8_32x32        , "cmatcholbkwsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholbkwsubst10x10_32x32      , "cmatcholbkwsubst10x10_32x32_cond20dynRange29qRAm3_0.seq" );
        DO_TEST( &cmatcholmmsesolver4x4_32x32      , "cmatcholmmsesolver4x4_32x32_cond2dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholmmsesolver4x4_32x32      , "cmatcholmmsesolver4x4_32x32_cond20dynRange29qRAm3_0.seq"  );
        DO_TEST( &cmatcholmmsesolver6x6_32x32      , "cmatcholmmsesolver6x6_32x32_cond2dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholmmsesolver6x6_32x32      , "cmatcholmmsesolver6x6_32x32_cond20dynRange29qRAm3_0.seq"  );
        DO_TEST( &cmatcholmmsesolver8x8_32x32      , "cmatcholmmsesolver8x8_32x32_cond2dynRange29qRAm3_0.seq"   );
        DO_TEST( &cmatcholmmsesolver8x8_32x32      , "cmatcholmmsesolver8x8_32x32_cond20dynRange29qRAm3_0.seq"  );
        DO_TEST( &cmatcholmmsesolver10x10_32x32    , "cmatcholmmsesolver10x10_32x32_cond2dynRange29qRAm3_0.seq" );
        DO_TEST( &cmatcholmmsesolver10x10_32x32    , "cmatcholmmsesolver10x10_32x32_cond20dynRange29qRAm3_0.seq");
        DO_TEST( &cmatcholpreprocess4x4_32x32      , "cmatcholpreprocess4x4_32x32_cond20dynRange29qRAm3_0.seq"      );
        DO_TEST( &cmatcholpreprocess6x6_32x32      , "cmatcholpreprocess6x6_32x32_cond20dynRange29qRAm3_0.seq"      );
        DO_TEST( &cmatcholpreprocess8x8_32x32      , "cmatcholpreprocess8x8_32x32_cond20dynRange29qRAm3_0.seq"      );
        DO_TEST( &cmatcholpreprocess10x10_32x32    , "cmatcholpreprocess10x10_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST(&cmatcholpseudoinv4x4_32x32        , "cmatcholpseudoinv4x4_32x32_cond2.seq"     );
        DO_TEST(&cmatcholpseudoinv6x6_32x32        , "cmatcholpseudoinv6x6_32x32_cond2.seq"     );
        DO_TEST(&cmatcholpseudoinv8x8_32x32        , "cmatcholpseudoinv8x8_32x32_cond2.seq"     );
        DO_TEST(&cmatcholpseudoinv10x10_32x32      , "cmatcholpseudoinv10x10_32x32_cond2.seq"   );
        DO_TEST(&cmatcholpseudoinv4x4_32x32        , "cmatcholpseudoinv4x4_32x32_cond20.seq"    );
        DO_TEST(&cmatcholpseudoinv6x6_32x32        , "cmatcholpseudoinv6x6_32x32_cond20.seq"    );
        DO_TEST(&cmatcholpseudoinv8x8_32x32        , "cmatcholpseudoinv8x8_32x32_cond20.seq"    );
        DO_TEST(&cmatcholpseudoinv10x10_32x32      , "cmatcholpseudoinv10x10_32x32_cond20.seq"  );
        DO_TEST( &matcholdecomp4x4_32x32           , "matcholdecomp4x4_32x32_cond20dynRange29qRAm3_0.seq"       );
        DO_TEST( &matcholdecomp4x4_32x32           , "matcholdecomp4x4_32x32_cond200dynRange29qRAm3_0.seq"      );
        DO_TEST( &matcholdecomp4x4_32x32           , "matcholdecomp4x4_32x32_cond2000dynRange29qRAm3_0.seq"     );
        DO_TEST( &matcholdecomp6x6_32x32           , "matcholdecomp6x6_32x32_cond20dynRange29qRAm3_0.seq"       );
        DO_TEST( &matcholdecomp6x6_32x32           , "matcholdecomp6x6_32x32_cond200dynRange29qRAm3_0.seq"      );
        DO_TEST( &matcholdecomp6x6_32x32           , "matcholdecomp6x6_32x32_cond2000dynRange29qRAm3_0.seq"     );
        DO_TEST( &matcholdecomp8x8_32x32           , "matcholdecomp8x8_32x32_cond20dynRange29qRAm3_0.seq"       );
        DO_TEST( &matcholdecomp8x8_32x32           , "matcholdecomp8x8_32x32_cond200dynRange29qRAm3_0.seq"      );
        DO_TEST( &matcholdecomp8x8_32x32           , "matcholdecomp8x8_32x32_cond2000dynRange29qRAm3_0.seq"     );
        DO_TEST( &matcholdecomp10x10_32x32         , "matcholdecomp10x10_32x32_cond20dynRange29qRAm3_0.seq"     );
        DO_TEST( &matcholdecomp10x10_32x32         , "matcholdecomp10x10_32x32_cond200dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholdecomp10x10_32x32         , "matcholdecomp10x10_32x32_cond2000dynRange29qRAm3_0.seq"   );
        DO_TEST( &matcholfwdsubst4x4_32x32         , "matcholfwdsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholfwdsubst6x6_32x32         , "matcholfwdsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholfwdsubst8x8_32x32         , "matcholfwdsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholfwdsubst10x10_32x32       , "matcholfwdsubst10x10_32x32_cond20dynRange29qRAm3_0.seq"  );
        DO_TEST( &matcholbkwsubst4x4_32x32         , "matcholbkwsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholbkwsubst6x6_32x32         , "matcholbkwsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholbkwsubst8x8_32x32         , "matcholbkwsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"    );
        DO_TEST( &matcholbkwsubst10x10_32x32       , "matcholbkwsubst10x10_32x32_cond20dynRange29qRAm3_0.seq"  );
        DO_TEST(&matcholmmsesolver4x4_32x32        , "matcholmmsesolver4x4_32x32_cond2dynRange29qRAm3_0.seq"    );
        DO_TEST(&matcholmmsesolver4x4_32x32        , "matcholmmsesolver4x4_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST(&matcholmmsesolver6x6_32x32        , "matcholmmsesolver6x6_32x32_cond2dynRange29qRAm3_0.seq"    );
        DO_TEST(&matcholmmsesolver6x6_32x32        , "matcholmmsesolver6x6_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST(&matcholmmsesolver8x8_32x32        , "matcholmmsesolver8x8_32x32_cond2dynRange29qRAm3_0.seq"    );
        DO_TEST(&matcholmmsesolver8x8_32x32        , "matcholmmsesolver8x8_32x32_cond20dynRange29qRAm3_0.seq"   );
        DO_TEST(&matcholmmsesolver10x10_32x32      , "matcholmmsesolver10x10_32x32_cond2dynRange29qRAm3_0.seq"  );
        DO_TEST(&matcholmmsesolver10x10_32x32      , "matcholmmsesolver10x10_32x32_cond20dynRange29qRAm3_0.seq" );
        DO_TEST( &matcholpreprocess4x4_32x32       , "matcholpreprocess4x4_32x32_cond20dynRange29qRAm3_0.seq"       );
        DO_TEST( &matcholpreprocess6x6_32x32       , "matcholpreprocess6x6_32x32_cond20dynRange29qRAm3_0.seq"       );
        DO_TEST( &matcholpreprocess8x8_32x32       , "matcholpreprocess8x8_32x32_cond20dynRange29qRAm3_0.seq"       );
        DO_TEST( &matcholpreprocess10x10_32x32     , "matcholpreprocess10x10_32x32_cond20dynRange29qRAm3_0.seq"     );
        DO_TEST(&matcholpseudoinv4x4_32x32         , "matcholpseudoinv4x4_32x32_cond2.seq"     );
        DO_TEST(&matcholpseudoinv6x6_32x32         , "matcholpseudoinv6x6_32x32_cond2.seq"     );
        DO_TEST(&matcholpseudoinv8x8_32x32         , "matcholpseudoinv8x8_32x32_cond2.seq"     );
        DO_TEST(&matcholpseudoinv10x10_32x32       , "matcholpseudoinv10x10_32x32_cond2.seq"   );
        DO_TEST(&matcholpseudoinv4x4_32x32         , "matcholpseudoinv4x4_32x32_cond20.seq"     );
        DO_TEST(&matcholpseudoinv6x6_32x32         , "matcholpseudoinv6x6_32x32_cond20.seq"     );
        DO_TEST(&matcholpseudoinv8x8_32x32         , "matcholpseudoinv8x8_32x32_cond20.seq"     );
        DO_TEST(&matcholpseudoinv10x10_32x32       , "matcholpseudoinv10x10_32x32_cond20.seq"   );

  }
  if (phaseNum==0 || phaseNum==2)
  {
        DO_TEST( &cmatcholdecomp4x4f        , "cmatcholdecomp4x4nf_cond2.seq"          );
        DO_TEST( &cmatcholdecomp4x4f        , "cmatcholdecomp4x4nf_cond10.seq"         );
        DO_TEST( &cmatcholdecomp6x6f        , "cmatcholdecomp6x6nf_cond2.seq"          );
        DO_TEST( &cmatcholdecomp6x6f        , "cmatcholdecomp6x6nf_cond10.seq"         );
        DO_TEST( &cmatcholdecomp8x8f        , "cmatcholdecomp8x8nf_cond2.seq"          );
        DO_TEST( &cmatcholdecomp8x8f        , "cmatcholdecomp8x8nf_cond10.seq"         );
        DO_TEST( &cmatcholdecomp10x10f      , "cmatcholdecomp10x10nf_cond2.seq"        );
        DO_TEST( &cmatcholdecomp10x10f      , "cmatcholdecomp10x10nf_cond10.seq"       );

        DO_TEST( &cmatcholfwdsubst4x4f      , "cmatcholfwdsubst4x4x1nf_cond2.seq"          );
        DO_TEST( &cmatcholfwdsubst4x4f      , "cmatcholfwdsubst4x4x1nf_cond10.seq"         );
        DO_TEST( &cmatcholfwdsubst6x6f      , "cmatcholfwdsubst6x6x1nf_cond2.seq"          );
        DO_TEST( &cmatcholfwdsubst6x6f      , "cmatcholfwdsubst6x6x1nf_cond10.seq"         );
        DO_TEST( &cmatcholfwdsubst8x8f      , "cmatcholfwdsubst8x8x1nf_cond2.seq"          );
        DO_TEST( &cmatcholfwdsubst8x8f      , "cmatcholfwdsubst8x8x1nf_cond10.seq"         );
        DO_TEST( &cmatcholfwdsubst10x10f    , "cmatcholfwdsubst10x10x1nf_cond2.seq"        );
        DO_TEST( &cmatcholfwdsubst10x10f    , "cmatcholfwdsubst10x10x1nf_cond10.seq"       );

        DO_TEST( &cmatcholbkwsubst4x4f      , "cmatcholbkwsubst4x1nf_cond2.seq"          );
        DO_TEST( &cmatcholbkwsubst4x4f      , "cmatcholbkwsubst4x1nf_cond10.seq"         );
        DO_TEST( &cmatcholbkwsubst6x6f      , "cmatcholbkwsubst6x1nf_cond2.seq"          );
        DO_TEST( &cmatcholbkwsubst6x6f      , "cmatcholbkwsubst6x1nf_cond10.seq"         );
        DO_TEST( &cmatcholbkwsubst8x8f      , "cmatcholbkwsubst8x1nf_cond2.seq"          );
        DO_TEST( &cmatcholbkwsubst8x8f      , "cmatcholbkwsubst8x1nf_cond10.seq"         );
        DO_TEST( &cmatcholbkwsubst10x10f    , "cmatcholbkwsubst10x1nf_cond2.seq"         );
        DO_TEST( &cmatcholbkwsubst10x10f    , "cmatcholbkwsubst10x1nf_cond10.seq"        );

        DO_TEST(&cmatcholmmsesolver4x4f      , "cmatcholmmsesolver4x4x1nf_cond2.seq"   );
        DO_TEST(&cmatcholmmsesolver4x4f      , "cmatcholmmsesolver4x4x1nf_cond10.seq"  );
        DO_TEST(&cmatcholmmsesolver6x6f      , "cmatcholmmsesolver6x6x1nf_cond2.seq"   );
        DO_TEST(&cmatcholmmsesolver6x6f      , "cmatcholmmsesolver6x6x1nf_cond10.seq"  );
        DO_TEST(&cmatcholmmsesolver8x8f      , "cmatcholmmsesolver8x8x1nf_cond2.seq"   );
        DO_TEST(&cmatcholmmsesolver8x8f      , "cmatcholmmsesolver8x8x1nf_cond10.seq"  );
        DO_TEST(&cmatcholmmsesolver10x10f    , "cmatcholmmsesolver10x10x1nf_cond2.seq" );
        DO_TEST(&cmatcholmmsesolver10x10f    , "cmatcholmmsesolver10x10x1nf_cond10.seq");

        DO_TEST(&cmatcholpreprocess4x4f      , "cmatcholpreprocess4x4nf_cond2.seq"   );
        DO_TEST(&cmatcholpreprocess4x4f      , "cmatcholpreprocess4x4nf_cond10.seq"  );
        DO_TEST(&cmatcholpreprocess6x6f      , "cmatcholpreprocess6x6nf_cond2.seq"   );
        DO_TEST(&cmatcholpreprocess6x6f      , "cmatcholpreprocess6x6nf_cond10.seq"  );
        DO_TEST(&cmatcholpreprocess8x8f      , "cmatcholpreprocess8x8nf_cond2.seq"   );
        DO_TEST(&cmatcholpreprocess8x8f      , "cmatcholpreprocess8x8nf_cond10.seq"  );
        DO_TEST(&cmatcholpreprocess10x10f    , "cmatcholpreprocess10x10nf_cond2.seq" );
        DO_TEST(&cmatcholpreprocess10x10f    , "cmatcholpreprocess10x10nf_cond10.seq");

        DO_TEST(&cmatcholpseudoinv4x4f      , "cmatcholpseudoinv4x4nf_cond2.seq"   );
        DO_TEST(&cmatcholpseudoinv4x4f      , "cmatcholpseudoinv4x4nf_cond10.seq"  );
        DO_TEST(&cmatcholpseudoinv6x6f      , "cmatcholpseudoinv6x6nf_cond2.seq"   );
        DO_TEST(&cmatcholpseudoinv6x6f      , "cmatcholpseudoinv6x6nf_cond10.seq"  );
        DO_TEST(&cmatcholpseudoinv8x8f      , "cmatcholpseudoinv8x8nf_cond2.seq"   );
        DO_TEST(&cmatcholpseudoinv8x8f      , "cmatcholpseudoinv8x8nf_cond10.seq"  );
        DO_TEST(&cmatcholpseudoinv10x10f    , "cmatcholpseudoinv10x10nf_cond2.seq" );
        DO_TEST(&cmatcholpseudoinv10x10f    , "cmatcholpseudoinv10x10nf_cond10.seq");

        DO_TEST( &matcholdecomp4x4f      , "matcholdecomp4x4nf_cond2.seq"          );
        DO_TEST( &matcholdecomp4x4f      , "matcholdecomp4x4nf_cond10.seq"         );
        DO_TEST( &matcholdecomp6x6f      , "matcholdecomp6x6nf_cond2.seq"          );
        DO_TEST( &matcholdecomp6x6f      , "matcholdecomp6x6nf_cond10.seq"         );
        DO_TEST( &matcholdecomp8x8f      , "matcholdecomp8x8nf_cond2.seq"          );
        DO_TEST( &matcholdecomp8x8f      , "matcholdecomp8x8nf_cond10.seq"         );
        DO_TEST( &matcholdecomp10x10f    , "matcholdecomp10x10nf_cond2.seq"        );
        DO_TEST( &matcholdecomp10x10f    , "matcholdecomp10x10nf_cond10.seq"       );

        DO_TEST( &matcholfwdsubst4x4f      , "matcholfwdsubst4x4x1nf_cond2.seq"   );
        DO_TEST( &matcholfwdsubst4x4f      , "matcholfwdsubst4x4x1nf_cond10.seq"  );
        DO_TEST( &matcholfwdsubst6x6f      , "matcholfwdsubst6x6x1nf_cond2.seq"   );
        DO_TEST( &matcholfwdsubst6x6f      , "matcholfwdsubst6x6x1nf_cond10.seq"  );
        DO_TEST( &matcholfwdsubst8x8f      , "matcholfwdsubst8x8x1nf_cond2.seq"   );
        DO_TEST( &matcholfwdsubst8x8f      , "matcholfwdsubst8x8x1nf_cond10.seq"  );
        DO_TEST( &matcholfwdsubst10x10f    , "matcholfwdsubst10x10x1nf_cond2.seq" );
        DO_TEST( &matcholfwdsubst10x10f    , "matcholfwdsubst10x10x1nf_cond10.seq");

        DO_TEST( &matcholbkwsubst4x4f      , "matcholbkwsubst4x1nf_cond2.seq"   );
        DO_TEST( &matcholbkwsubst4x4f      , "matcholbkwsubst4x1nf_cond10.seq"  );
        DO_TEST( &matcholbkwsubst6x6f      , "matcholbkwsubst6x1nf_cond2.seq"   );
        DO_TEST( &matcholbkwsubst6x6f      , "matcholbkwsubst6x1nf_cond10.seq"  );
        DO_TEST( &matcholbkwsubst8x8f      , "matcholbkwsubst8x1nf_cond2.seq"   );
        DO_TEST( &matcholbkwsubst8x8f      , "matcholbkwsubst8x1nf_cond10.seq"  );
        DO_TEST( &matcholbkwsubst10x10f    , "matcholbkwsubst10x1nf_cond2.seq"  );
        DO_TEST( &matcholbkwsubst10x10f    , "matcholbkwsubst10x1nf_cond10.seq" );
	  
        DO_TEST(&matcholmmsesolver4x4f      , "matcholmmsesolver4x4x1nf_cond2.seq"   );
        DO_TEST(&matcholmmsesolver4x4f      , "matcholmmsesolver4x4x1nf_cond10.seq"  );
        DO_TEST(&matcholmmsesolver6x6f      , "matcholmmsesolver6x6x1nf_cond2.seq"   );
        DO_TEST(&matcholmmsesolver6x6f      , "matcholmmsesolver6x6x1nf_cond10.seq"  );
        DO_TEST(&matcholmmsesolver8x8f      , "matcholmmsesolver8x8x1nf_cond2.seq"   );
        DO_TEST(&matcholmmsesolver8x8f      , "matcholmmsesolver8x8x1nf_cond10.seq"  );
        DO_TEST(&matcholmmsesolver10x10f    , "matcholmmsesolver10x10x1nf_cond2.seq" );
        DO_TEST(&matcholmmsesolver10x10f    , "matcholmmsesolver10x10x1nf_cond10.seq");

        DO_TEST(&matcholpreprocess4x4f      , "matcholpreprocess4x4nf_cond2.seq"   );
        DO_TEST(&matcholpreprocess4x4f      , "matcholpreprocess4x4nf_cond10.seq"  );
        DO_TEST(&matcholpreprocess6x6f      , "matcholpreprocess6x6nf_cond2.seq"   );
        DO_TEST(&matcholpreprocess6x6f      , "matcholpreprocess6x6nf_cond10.seq"  );
        DO_TEST(&matcholpreprocess8x8f      , "matcholpreprocess8x8nf_cond2.seq"   );
        DO_TEST(&matcholpreprocess8x8f      , "matcholpreprocess8x8nf_cond10.seq"  );
        DO_TEST(&matcholpreprocess10x10f    , "matcholpreprocess10x10nf_cond2.seq" );
        DO_TEST(&matcholpreprocess10x10f    , "matcholpreprocess10x10nf_cond10.seq");

        DO_TEST(&matcholpseudoinv4x4f  , "matcholpseudoinv4x4nf_cond2.seq");
        DO_TEST(&matcholpseudoinv4x4f  , "matcholpseudoinv4x4nf_cond10.seq");
        DO_TEST(&matcholpseudoinv6x6f  , "matcholpseudoinv6x6nf_cond2.seq");
        DO_TEST(&matcholpseudoinv6x6f  , "matcholpseudoinv6x6nf_cond10.seq");
        DO_TEST(&matcholpseudoinv8x8f  , "matcholpseudoinv8x8nf_cond2.seq");
        DO_TEST(&matcholpseudoinv8x8f  , "matcholpseudoinv8x8nf_cond10.seq");
        DO_TEST(&matcholpseudoinv10x10f, "matcholpseudoinv10x10nf_cond2.seq");
        DO_TEST(&matcholpseudoinv10x10f, "matcholpseudoinv10x10nf_cond10.seq");

  }

  return res;
}
