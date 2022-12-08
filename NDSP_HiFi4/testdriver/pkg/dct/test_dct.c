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
 * Test procedures for real DCT functions.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(fft)
/* Test engine API. */
#include "testeng.h"
/* Test engine extension for DCT. */
#include "testeng_dct.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"

/* Suppress Visual C warnings on +/-HUGE_VALF macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

#define sz_fr16c    sizeof(complex_fract16)
#define sz_fl64c    sizeof(complex_double)

/* Initializer for a test description structure. */
#define TEST_DESC_STND( fmt, extraParam, align,fwddct,invdct )   {{ (fmt),(extraParam),NULL,TE_DIM_NUM_1,(align), \
                                                     &te_create_dct, &te_destroy_dct, \
                                                     &te_load_dct, &fileProc_stnd },(tFrwTransFxn)fwddct,(tInvTransFxn)invdct}
#define TEST_DESC_DCT2D( fmt, extraParam, align,fwddct,invdct )   {{ (fmt),(extraParam),NULL,TE_DIM_NUM_2,(align), \
                                                     &te_create_dct, &te_destroy_dct, \
                                                     &te_loadFxn_dct2d, &te_processFxn_dct2d },(tFrwTransFxn)fwddct,(tInvTransFxn)invdct}



static const struct tagTestDef
{
  int                phaseNum;
  const char *       seqFilename;
  tTestEngTarget     target;
  tTestEngDesc_dct   desc;

} testTbl[] = 
{
  { 1, "dct32_16x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT16, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_16x16 , NULL   ) },
  { 1, "dct64_16x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT16, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_16x16 , NULL   ) },
#if 0 // HiFi3/3z API
  { 1, "dct32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_24x24 , NULL   ) },
  { 1, "dct64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_24x24 , NULL   ) },
#endif
  { 1, "dct32_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x16 , NULL   ) },
  { 1, "dct64_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x16 , NULL   ) },
  { 1, "dct32_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x32 , NULL   ) },
  { 1, "dct64_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x32 , NULL   ) },
#if 0 // HiFi3/3z API
  { 1, "dct4_32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_128_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_256_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_512_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_128_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_256_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct4_512_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
#endif
  { 1, "dct4_32_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct4_64_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct4_128_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct4_256_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct4_512_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },

  { 1, "dct4_32_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct4_64_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct4_128_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct4_256_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct4_512_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 2, "dct32f.seq"      , (tTestEngTarget)te_frameProc_stnd_dct     , TEST_DESC_STND( FMT_FLOAT32, ATTR_FAST                      , TE_ALIGN_YES, dctf      , NULL   ) },
  { 2, "dct64f.seq"      , (tTestEngTarget)te_frameProc_stnd_dct     , TEST_DESC_STND( FMT_FLOAT32, ATTR_FAST                      , TE_ALIGN_YES, dctf      , NULL   ) },

#if 0 // HiFi3/3z API
  { 1, "mdct32_24x24.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "mdct64_24x24.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "mdct128_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "mdct256_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "mdct512_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
#endif
  { 1, "mdct32_32x16.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "mdct64_32x16.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "mdct128_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "mdct256_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "mdct512_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },

  { 1, "mdct32_32x32.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "mdct64_32x32.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "mdct128_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "mdct256_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "mdct512_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
#if 0 // HiFi3/3z API
  { 1, "imdct32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "imdct64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "imdct128_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "imdct256_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "imdct512_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
#endif
  { 1, "imdct32_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "imdct64_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "imdct128_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "imdct256_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "imdct512_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },

  { 1, "imdct32_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "imdct64_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "imdct128_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "imdct256_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "imdct512_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },

  { 1, "dct2d_8x16.seq" , (tTestEngTarget)dct2d_8x16 , TEST_DESC_DCT2D( FMT_UINT8|FMT_REAL, TE_FFT_OPT_SCALE_METH|TE_DCT2D , TE_ALIGN_YES, dct2d_8x16, NULL ) },
  { 1, "idct2d_16x8.seq", (tTestEngTarget)idct2d_16x8, TEST_DESC_DCT2D( FMT_INT16|FMT_REAL, TE_FFT_OPT_SCALE_METH|TE_IDCT2D, TE_ALIGN_YES, NULL, idct2d_16x8) },

  { 0 } /* End of table */
};
 
/* Perform all tests for dct API functions. */
int main_dct( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;

  for ( ix=0,res=1; testTbl[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl[ix].phaseNum )
    {
      res &= ( 0 != TestEngRun( (tTestEngTarget )testTbl[ix].target,
                                &testTbl[ix].desc.desc,
                                testTbl[ix].seqFilename,
                                isFull, isVerbose, breakOnError,0 ) );
    }
  }

  return (res);

} /* main_dct() */
