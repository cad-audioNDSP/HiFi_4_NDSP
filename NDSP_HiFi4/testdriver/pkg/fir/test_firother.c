/* ------------------------------------------------------------------------ */
/* Copyright (c) 2018 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ("Cadence    */
/* Libraries") are the copyrighted works of Cadence Design Systems Inc.	    */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
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
/*          Copyright (C) 2015-2018 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
 * Test procedures for FIR
 */

#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API */
#include LIBRARY_HEADER(fir)
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
#include "testeng_fir.h"

#define MAX_FUNC_NUM  4
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extraParam, dimNum, align, loadFxn, procFxn ) { (fmt),extraParam,NULL, (dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }
#define TEST_DESC_CONVOLVE( fmt, extraParam, align )  { (fmt), (extraParam)|TE_FIR_CONVOLVE_API , NULL, TE_DIM_NUM_2, (align), NULL,NULL, &te_loadFxn_crosscorr, &te_processFxn_crosscorr }
#define TEST_DESC_CROSSCORR( fmt, extraParam, align ) { (fmt), (extraParam)|TE_FIR_CROSSCORR_API, NULL, TE_DIM_NUM_2, (align), NULL,NULL, &te_loadFxn_crosscorr, &te_processFxn_crosscorr }
#define TEST_DESC_AUTOCORR( fmt, extraParam, align )  { (fmt), (extraParam)|TE_FIR_AUTOCORR_API , NULL, TE_DIM_NUM_1, (align), NULL,NULL, &te_loadFxn_autocorr,  &te_processFxn_autocorr  }

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/*-----------------------------------------------
wrapper for LMS to test convergence:
input:
x[M+N*P-1]  far end
r[N*P]      near end
norm[P]     norm factor for P blocks
mu          adapt rate
M           IR length
N           block size
P           number of blocks
output:
h[M]        estimated impulse response
temporary
e[N]        error
-----------------------------------------------*/
void fir_blmsf_convergence( 
                      float32_t * e, float32_t * h, const float32_t * r,
                const float32_t * x, const float32_t* norm, float32_t mu, int N, int M, int P )
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blmsf(e,h,r,x,*norm,mu,N,M);
    }
}

void cxfir_blmsf_convergence( complex_float * e, complex_float * h, 
                const complex_float * r,
                const complex_float * x, 
                const float32_t* norm, float32_t mu, 
                int          N, int       M , int P)
{
    int m,p;
    for(m=0; m<2*M; m++) ((float32_t*)h)[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        cxfir_blmsf(e,h,r,x,*norm,mu,N,M);
    }
}
void fir_blms16x16_convergence (  int16_t * e, int16_t *  h,
                const int16_t * r,
                const int16_t * x,
                const int16_t * norm,int16_t   mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms16x16(e,h,r,x,*norm,mu,N,M);
    }
}
void fir_blms16x32_convergence (  int32_t * e, int32_t *  h,
                const int16_t * r,
                const int16_t * x,
                const int32_t * norm,int16_t   mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms16x32(e,h,r,x,*norm,mu,N,M);
    }
}
void fir_blms32x32_convergence(  int32_t * e, int32_t *  h,
                const int32_t * r,
                const int32_t * x,
                const int32_t * norm, int32_t mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms32x32(e,h,r,x,*norm,mu,N,M);
    }
}
void fir_blms32x32ep_convergence(  int32_t * e, int32_t *  h,
                const int32_t * r,
                const int32_t * x,
                const int32_t * norm, int32_t mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms32x32ep(e,h,r,x,*norm,mu,N,M);
    }
}
void cxfir_blms32x32_convergence (complex_fract32 *  e, complex_fract32 *  h,
                const complex_fract32 *  r,
                const complex_fract32 *  x,
                const int32_t *  norm, int32_t mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<2*M; m++) ((int32_t*)h)[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        cxfir_blms32x32(e,h,r,x,*norm,mu,N,M);
    }
}

/* check if target function is not available in the target configuration already */
static int te_create_convergence(tTestEngContext * context)
{
    return IS_PRESENT(context->desc->extraPtr) ? 1:-1;
}

  /* API test definitions. */
  static const struct 
  {
    tTestEngTarget   funcList[MAX_FUNC_NUM];
    tTestEngDesc     testDesc;
  }
  testDefTbl[] =
  {
    /* Stage 1 */
    {  FUNC_LIST( (tTestEngTarget)&fir_convol32x32,
                  (tTestEngTarget)&fir_convol32x32ep),  TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convol16x16),    TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_YES ) },
#if 0// for HiFi3/3z
    {  FUNC_LIST( (tTestEngTarget)&fir_convol24x24),    TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_24X24, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola24x24),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_24X24, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorr24x24),     TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_24X24, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorra24x24),    TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_24X24, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr24x24),     TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_24X24, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra24x24),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_24X24, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms24x24),      { FMT_REAL|FMT_FRACT32, 0                 ,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
#endif
    {  FUNC_LIST( (tTestEngTarget)&fir_convol32x16),    TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_convol32x16),  TEST_DESC_CONVOLVE( FMT_CPLX|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola16x16),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola32x32),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola32x32ep), TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_EP   , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola32x16),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_convola32x16), TEST_DESC_CONVOLVE( FMT_CPLX|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lconvola16x16),  TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lconvola32x32),  TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorr32x32,
                  (tTestEngTarget)&fir_acorr32x32ep),   TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorr16x16),     TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorra16x16),    TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorra32x32),    TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorra32x32ep),  TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_EP   , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lacorra16x16),   TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT16, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lacorra32x32),   TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr32x32,
                  (tTestEngTarget)&fir_xcorr32x32ep),   TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorr32x32),   TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr16x16),     TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr32x16),     TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra16x16),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra32x32),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorra32x32),  TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra32x32ep),  TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, TE_FIR_OTHER_EP   , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra32x16),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lxcorra16x16),   TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lxcorra32x32),   TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },

    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32,
                  (tTestEngTarget)&fir_blms32x32ep),    { FMT_REAL|FMT_FRACT32, 0                 ,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x16),      { FMT_REAL|FMT_FRACT16, 0                 ,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x32),      { FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blms32x32),    { FMT_CPLX|FMT_FRACT32, 0                 ,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },

    /* Stage 2 */
    {  FUNC_LIST( (tTestEngTarget)&fir_convolf),        TEST_DESC_CONVOLVE ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorrf),         TEST_DESC_AUTOCORR ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorrf),         TEST_DESC_CROSSCORR( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorrf),       TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convolaf),       TEST_DESC_CONVOLVE ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorraf),        TEST_DESC_AUTOCORR ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorraf),        TEST_DESC_CROSSCORR( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorraf),      TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_blmsf),          { FMT_REAL|FMT_FLOAT32, 0,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blmsf),        { FMT_CPLX|FMT_FLOAT32, 0,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },

    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x16_convergence),      { FMT_REAL|FMT_FRACT16, 0                 ,(void *)fir_blms16x16  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x32_convergence),      { FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16,(void *)fir_blms16x32  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32_convergence),      { FMT_REAL|FMT_FRACT32, 0,                 (void *)fir_blms32x32  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32ep_convergence),    { FMT_REAL|FMT_FRACT32, 0,                 (void *)fir_blms32x32ep,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blms32x32_convergence),    { FMT_CPLX|FMT_FRACT32, 0,                 (void *)cxfir_blms32x32,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blmsf_convergence),          { FMT_REAL|FMT_FLOAT32, 0,                 (void *)fir_blmsf      ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blmsf_convergence),        { FMT_CPLX|FMT_FLOAT32, 0,                 (void *)cxfir_blmsf    ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },

    {  FUNC_LIST( NULL ),              TEST_DESC(  0, 0, 0, 0, NULL, NULL ) } /* End of table */
  };



/* Perform all tests for FIR API functions. */
int main_firother( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
 
    int res = 1;

    #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), (seqFile),         \
                                                       isFull, isVerbose, breakOnError ) )
    if ( phaseNum == 0 || phaseNum == 1 )
    {
        DO_TEST( &fir_acorr16x16     , "fir_acorr16x16.seq"      );
#if 0// for HiFi3/3z
        DO_TEST( &fir_acorr24x24     , "fir_acorr24x24.seq"      );
#endif
        DO_TEST( &fir_acorr32x32     , "fir_acorr32x32.seq"      );
        DO_TEST( &fir_acorr32x32ep   , "fir_acorr32x32ep.seq"    );
        DO_TEST( &fir_xcorr16x16     , "fir_xcorr16x16.seq"      );
        DO_TEST( &fir_xcorr32x16     , "fir_xcorr32x16.seq"      );
#if 0// for HiFi3/3z
        DO_TEST( &fir_xcorr24x24     , "fir_xcorr24x24.seq"      );
#endif
        DO_TEST( &fir_xcorr32x32     , "fir_xcorr32x32.seq"      );
        DO_TEST( &fir_xcorr32x32ep   , "fir_xcorr32x32ep.seq"    );
        DO_TEST( &cxfir_xcorr32x32   , "cxfir_xcorr32x32.seq"    );
        DO_TEST( &cxfir_xcorra32x32  , "cxfir_xcorra32x32.seq"     );
        DO_TEST( &fir_convol16x16    , "fir_convol16x16.seq"     );
        DO_TEST( &fir_convol32x16    , "fir_convol32x16.seq"     );
#if 0// for HiFi3/3z
        DO_TEST( &fir_convol24x24    , "fir_convol24x24.seq"     );
#endif
        DO_TEST( &fir_convol32x32    , "fir_convol32x32.seq"     );
        DO_TEST( &fir_convol32x32ep  , "fir_convol32x32ep.seq"   );
        DO_TEST( &cxfir_convol32x16  , "cxfir_convol32x16.seq"   );
        DO_TEST( &fir_acorra16x16    , "fir_acorra16x16.seq"     );
#if 0// for HiFi3/3z
        DO_TEST( &fir_acorra24x24    , "fir_acorra24x24.seq"     );
#endif
        DO_TEST( &fir_acorra32x32    , "fir_acorra32x32.seq"     );
        DO_TEST( &fir_acorra32x32ep  , "fir_acorra32x32ep.seq"   );
        DO_TEST( &fir_lacorra16x16   , "fir_lacorra16x16.seq"    );
        DO_TEST( &fir_lacorra32x32   , "fir_lacorra32x32.seq"    );
        DO_TEST( &fir_xcorra16x16    , "fir_xcorra16x16.seq"     );
        DO_TEST( &fir_xcorra32x16    , "fir_xcorra32x16.seq"     );
#if 0// for HiFi3/3z
        DO_TEST( &fir_xcorra24x24    , "fir_xcorra24x24.seq"     );
#endif
        DO_TEST( &fir_xcorra32x32    , "fir_xcorra32x32.seq"     );
        DO_TEST( &fir_xcorra32x32ep  , "fir_xcorra32x32ep.seq"   );
        DO_TEST( &fir_lxcorra16x16   , "fir_lxcorra16x16.seq"    );
        DO_TEST( &fir_lxcorra32x32   , "fir_lxcorra32x32.seq"    );
        DO_TEST( &fir_convola16x16   , "fir_convola16x16.seq"    );
        DO_TEST( &fir_convola32x16   , "fir_convola32x16.seq"    );
#if 0// for HiFi3/3z
        DO_TEST( &fir_convola24x24   , "fir_convola24x24.seq"    );
#endif
        DO_TEST( &fir_convola32x32   , "fir_convola32x32.seq"    );
        DO_TEST( &fir_convola32x32ep , "fir_convola32x32ep.seq"  );
        DO_TEST( &cxfir_convola32x16 , "cxfir_convola32x16.seq"  );
        DO_TEST( &fir_lconvola16x16  , "fir_lconvola16x16.seq"   );
        DO_TEST( &fir_lconvola32x32  , "fir_lconvola32x32.seq"   );

        DO_TEST( &fir_blms16x16      , "fir_blms16x16.seq"       );
#if 0// for HiFi3/3z
        DO_TEST( &fir_blms24x24      , "fir_blms24x24.seq"       );
#endif
        DO_TEST( &fir_blms32x32      , "fir_blms32x32.seq"       );
        DO_TEST( &fir_blms32x32ep    , "fir_blms32x32ep.seq"     );
        DO_TEST( &fir_blms16x32      , "fir_blms16x32.seq"       );
        DO_TEST( &cxfir_blms32x32    , "cxfir_blms32x32.seq"     );
        DO_TEST( &fir_blms16x16_convergence    , "fir_blms16x16_convergence.seq"    );
        DO_TEST( &fir_blms16x32_convergence    , "fir_blms16x32_convergence.seq"    );
        DO_TEST( &cxfir_blms32x32_convergence  , "cxfir_blms32x32_convergence.seq"  );
        DO_TEST( &fir_blms32x32_convergence    , "fir_blms32x32_convergence.seq"    );
        DO_TEST( &fir_blms32x32ep_convergence  , "fir_blms32x32ep_convergence.seq"  );
    }
    if ( phaseNum == 0 || phaseNum == 2)
    {
        DO_TEST( &fir_acorrf       , "fir_acorrf.seq"       );
        DO_TEST( &fir_xcorrf       , "fir_xcorrf.seq"       );
        DO_TEST( &cxfir_xcorrf     , "cxfir_xcorrf.seq"     );
        DO_TEST( &fir_convolf      , "fir_convolf.seq"      );
        DO_TEST( &fir_acorraf      , "fir_acorraf.seq"      );
        DO_TEST( &fir_xcorraf      , "fir_xcorraf.seq"      );
        DO_TEST( &cxfir_xcorraf    , "cxfir_xcorraf.seq"    );
        DO_TEST( &fir_convolaf     , "fir_convolaf.seq"     );
        DO_TEST( &fir_blmsf        , "fir_blmsf.seq"        );
        DO_TEST( &cxfir_blmsf      , "cxfir_blmsf.seq"      );
        DO_TEST( &cxfir_blmsf_convergence      , "cxfir_blmsf_convergence.seq"      );
        DO_TEST( &fir_blmsf_convergence        , "fir_blmsf_convergence.seq"        );
    }

    return (res);
} /* main_firother() */
