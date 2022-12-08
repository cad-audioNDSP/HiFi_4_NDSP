/* ------------------------------------------------------------------------ */
/* Copyright (c) 2018 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs (“Cadence    */
/* Libraries”) are the copyrighted works of Cadence Design Systems Inc.	    */
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
 * Test procedures for real FFT functions.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(fft)
/* Fixed point arithmetics. */
#include "NatureDSP_Math.h"
/* Test engine API. */
#include "testeng.h"
/* Test engine extension for FFT. */
#include "testeng_fft.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Suppress Visual C warnings on +/-HUGE_VALF macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

#define sz_fr16c    sizeof(complex_int16_t)
#define sz_fl64c    sizeof(complex_double)

/* Period, in frames, between reallocation of in/out buffers consumed by the
 * target FFT routine. */
#define XY_VEC_REALLOC_PERIOD     16

/* Test data file processing for a standalone test. Applies the target FFT function
 * to test data, compares the FFT output with reference data and estimates the SINAD. */
static void fileProc_stnd( tTestEngContext * context );

#ifdef PSEUDOFLOAT
//Forward FFT 
void float32_fft(int16_t *in, float64_t *out, int N);
void float64_fft(int16_t *in, float64_t *out, int N);
//Inverse FFT 
//void float32_ifft(int16_t* in, float64_t* out, int N);
//void float64_ifft(int16_t* in, float64_t* out, int N);
#endif

/* Initializer for a test description structure. */
#define TEST_DESC_STND( fmt, extraParam, align,fwdfft,invfft )   {{ (fmt),(extraParam),NULL, TE_ARGNUM_1,(align), \
                                                     &te_create_fft, &te_destroy_fft, \
                                                     &te_load_fft, &fileProc_stnd },(tFrwTransFxn)fwdfft,(tInvTransFxn)invfft}

/* FFT attributes for feature rich and fast real FFTs. */
#define ATTR_FEATURE_RICH      (TE_FFT_REAL|TE_FFT_OPT_INPLACE|TE_FFT_FULL_SPECTRUM)
#define ATTR_FEATURE_RICH_SCL  (TE_FFT_REAL|TE_FFT_OPT_INPLACE|TE_FFT_FULL_SPECTRUM|TE_FFT_OPT_SCALE_METH)
#define ATTR_FAST              (TE_FFT_REAL|TE_FFT_OPT_REUSE_INBUF)

typedef const struct tagTestDef
{
  int                phaseNum;
  const char *       seqFilename;
  tTestEngTarget     target;
  tTestEngDesc_fft   desc;

}TestDef_t;

static TestDef_t testTbl_rfft[] =
{
#if 1
    { 1, "fft32_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft64_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft128_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft256_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft512_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft1024_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft2048_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft4096_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft8192_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "fft32_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft64_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft128_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft256_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft512_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft1024_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft2048_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft4096_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft8192_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "ifft32_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft64_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft128_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft256_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft512_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft1024_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft2048_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft4096_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft8192_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },

    { 1, "ifft32_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft64_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft128_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft256_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft512_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft1024_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft2048_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft4096_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft8192_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
#endif
#if 0 //HiFi3/3z API
    { 1, "fft_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "fft_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "fft_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "fft_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "ifft_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "ifft_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "ifft_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "ifft_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
#endif 
#if 1
    { 1, "fft32_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft64_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft128_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft256_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft512_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft1024_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft2048_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft4096_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft8192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "fft32_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft64_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft128_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft256_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft512_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft1024_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft2048_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft4096_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft8192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "ifft32_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft64_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft128_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft256_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft512_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft1024_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft2048_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft4096_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft8192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },

    { 1, "ifft32_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft64_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft128_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft256_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft512_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft1024_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft2048_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft4096_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft8192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
#endif

#if 1
    { 1, "fft32_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft64_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft128_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft256_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft512_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1024_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft2048_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft4096_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft8192_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "fft32_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft64_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft128_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft256_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft512_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1024_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft2048_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft4096_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft8192_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "ifft32_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft64_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft128_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft256_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft512_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1024_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft2048_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft4096_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },

    { 1, "ifft32_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft64_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft128_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft256_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft512_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1024_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft2048_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft4096_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft8192_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
#endif 
    { 0 } /* End of table */
}; //testTbl_rfft


static TestDef_t testTbl_rnfft[] =
{
#if 0
    // Mixed radix rfft
    { 1, "fft160_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft192_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft240_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft320_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft384_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft480_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "fft160_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft192_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft240_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft320_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft384_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "fft480_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "ifft160_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft192_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft240_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft320_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft384_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft480_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },

    { 1, "ifft160_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft192_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft240_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft320_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft384_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "ifft480_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
#endif
#if 0
    // Mixed radix rfft
    { 1, "fft160_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft240_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft320_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft384_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft480_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "fft160_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft240_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft320_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft384_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "fft480_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "ifft160_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft240_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft320_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft384_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft480_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },

    { 1, "ifft160_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft240_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft320_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft384_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "ifft480_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
#endif
#if 1

    // Mixed radix rnfft
    { 1, "fft12_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft24_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft30_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft36_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft48_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft60_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft72_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft90_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft96_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft108_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft120_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft144_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft160_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft180_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft192_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft216_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft240_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft288_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft300_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft320_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft324_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft360_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft384_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft432_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft480_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft540_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft576_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft720_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft768_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft960_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1152_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1440_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1536_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1920_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "fft12_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft24_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft30_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft36_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft48_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft60_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft72_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft90_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft96_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft108_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft120_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft144_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft160_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft180_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft192_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft216_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft240_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft288_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft300_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft320_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft324_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft360_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft384_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft432_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft480_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft540_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft576_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft720_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft768_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft960_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1152_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1440_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1536_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "fft1920_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    // Mixed radix rnifft
    { 1, "ifft12_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft24_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft30_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft36_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft48_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft60_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft72_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft90_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft96_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft108_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft120_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft144_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft160_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft180_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft192_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft216_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft240_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft288_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft300_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft320_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft324_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft360_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft384_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft432_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft480_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft540_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft576_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft720_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft768_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft960_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1152_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1440_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1536_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1920_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },

    { 1, "ifft12_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft24_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft30_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft36_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft48_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft60_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft72_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft90_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft96_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft108_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft120_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft144_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft160_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft180_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft192_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft216_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft240_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft288_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft300_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft320_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft324_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft360_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft384_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft432_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft480_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft540_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft576_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft720_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft768_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft960_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1152_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1440_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1536_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "ifft1920_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
#endif
    { 0 } /* End of table */
}; //testTbl_rnfft



static TestDef_t testTbl_rfftie[] =
{
#if 1
    { 1, "fft256_real16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, fft_real16x16_ie, NULL) },
    { 1, "fft512_real16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, fft_real16x16_ie, NULL) },
    { 1, "fft1024_real16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, fft_real16x16_ie, NULL) },
    { 1, "ifft256_real16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, ifft_real16x16_ie) },
    { 1, "ifft512_real16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, ifft_real16x16_ie) },
    { 1, "ifft1024_real16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, ifft_real16x16_ie) },
#endif
#if 0 //HiFi3/3z API
    { 1, "fft_real24x24_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, fft_real24x24_ie, NULL) },
    { 1, "ifft_real24x24_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_real24x24_ie) },
#endif
#if 1
    { 1, "fft256_real32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, fft_real32x16_ie, NULL) },
    { 1, "fft512_real32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, fft_real32x16_ie, NULL) },
    { 1, "fft1024_real32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, fft_real32x16_ie, NULL) },
    { 1, "fft256_real32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, fft_real32x16_ie, NULL) },
    { 1, "fft512_real32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, fft_real32x16_ie, NULL) },
    { 1, "fft1024_real32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, fft_real32x16_ie, NULL) },
    { 1, "ifft256_real32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, NULL, &ifft_real32x16_ie) },
    { 1, "ifft512_real32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, NULL, &ifft_real32x16_ie) },
    { 1, "ifft1024_real32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, NULL, &ifft_real32x16_ie) },
    { 1, "ifft256_real32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, NULL, &ifft_real32x16_ie) },
    { 1, "ifft512_real32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, NULL, &ifft_real32x16_ie) },
    { 1, "ifft1024_real32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16, TE_ALIGN_YES, NULL, &ifft_real32x16_ie) },
#endif
#if 1
    { 1, "fft256_real32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32_ie, NULL) },
    { 1, "fft512_real32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32_ie, NULL) },
    { 1, "fft1024_real32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32_ie, NULL) },
    { 1, "fft256_real32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32_ie, NULL) },
    { 1, "fft512_real32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32_ie, NULL) },
    { 1, "fft1024_real32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32_ie, NULL) },
    { 1, "ifft256_real32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32_ie) },
    { 1, "ifft512_real32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32_ie) },
    { 1, "ifft1024_real32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32_ie) },
    { 1, "ifft256_real32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32_ie) },
    { 1, "ifft512_real32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32_ie) },
    { 1, "ifft1024_real32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32_ie) },
#endif
#if 1
    { 2, "fft8_realf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft16_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft32_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft64_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft128_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft256_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft512_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft1024_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft2048_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "fft4096_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },

    { 2, "ifft8_realf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft16_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft32_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft64_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft128_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft256_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft512_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft1024_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft2048_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "ifft4096_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
#endif
    { 0 } /* End of table */
}; //testTbl_rfftie

/* Perform all tests for fft_realMxN, ifft_realMxN API functions. */
int main_rfft( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;
  for ( ix=0,res=1; testTbl_rfft[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl_rfft[ix].phaseNum )
    {
        tTestEngTarget target = testTbl_rfft[ix].target;
        /* Make sure that all functions is present */
        if (!IS_PRESENT(testTbl_rfft[ix].desc.frwTransFxn) &&
            !IS_PRESENT(testTbl_rfft[ix].desc.invTransFxn))
        {
            target = (tTestEngTarget)testTbl_rfft[ix].desc.frwTransFxn;
        }

        res &= ( 0 != TestEngRun(  target,
                                   &testTbl_rfft[ix].desc.desc,
                                   testTbl_rfft[ix].seqFilename,
                                   isFull, isVerbose, breakOnError,0 ) );
    }
  }
  return (res);
} /* main_rfft() */

/* Perform all tests for fft_realMxN, ifft_realMxN API functions. */
int main_rnfft( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;
  for ( ix=0,res=1; testTbl_rnfft[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl_rnfft[ix].phaseNum )
    {
        tTestEngTarget target = testTbl_rnfft[ix].target;
        /* Make sure that all functions is present */
        if (!IS_PRESENT(testTbl_rnfft[ix].desc.frwTransFxn) &&
            !IS_PRESENT(testTbl_rnfft[ix].desc.invTransFxn))
        {
            target = (tTestEngTarget)testTbl_rnfft[ix].desc.frwTransFxn;
        }

        res &= ( 0 != TestEngRun(  target,
                                   &testTbl_rnfft[ix].desc.desc,
                                   testTbl_rnfft[ix].seqFilename,
                                   isFull, isVerbose, breakOnError,0 ) );
    }
  }
  return (res);
} /* main_rnfft() */

/* Perform all tests for fft_realMxN_ie, ifft_realMxN_ie API functions. */
int main_rfftie(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_rfftie[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_rfftie[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_rfftie[ix].target;
            /* Make sure that all functions is present */
            if (!IS_PRESENT(testTbl_rfftie[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_rfftie[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_rfftie[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_rfftie[ix].desc.desc,
                testTbl_rfftie[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);
} /* main_rfftie() */

/* Test data file processing for a standalone test. Applies the target FFT function
 * to test data, compares the FFT output with reference data and estimates the SINAD. */
void fileProc_stnd( tTestEngContext * context )
{
  tTestEngContext_fft           * context_fft;
  tTestEngFrameProcFxn_fft_stnd * frameProcFxn;

  tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
  tVec xVec, yVec;            /* In/out vectors for the target FFT routine. */
  FILE *fIn = 0, *fRef = 0;
  int res = 0;
  int lenIn=0, lenOut=0;
  int lenFread = 0;  /* Number of samples readed from file */
  int fmtX=0, fmtY=0;
  int n, N;

  NASSERT( context && context->desc );
  NASSERT( context->target.fut && context->target.handle );

  context_fft = (tTestEngContext_fft*)context->target.handle;

  /* FFT size */
  N = context->args.N_DIM_;

  /* Select length and format (real/complex) for all vectors. */
  switch ( context->args.caseType )
  {
  case TE_FFT_TESTCASE_FORWARD:
    lenIn  = N;
    lenFread = N; 
    lenOut = ( ( context->desc->extraParam & TE_FFT_FULL_SPECTRUM ) ? N : N/2+1 );
    fmtX   = FMT_REAL;
    fmtY   = FMT_CPLX;
    break;
  case TE_FFT_TESTCASE_INVERSE:
    lenIn  = N/2+1;
    lenFread = lenIn - 1; 
    lenOut = N;
    fmtX   = FMT_CPLX;
    fmtY   = FMT_REAL;
    break;
  default:
    ASSERT( 0 );
  }

  /* Preset an invalid SINAD for the test result. */
  *vecGetElem_fl32( &context->dataSet.Z, 0 ) = -HUGE_VALF;

  if ( context->isVerbose )
  {
    if ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH )
    {
      printf( "scale_mtd %d ", context_fft->scale_method );
    }

    printf( "%-40s  ", context_fft->fRefName );
  }

  memset( &inVec , 0, sizeof(inVec ) );
  memset( &outVec, 0, sizeof(outVec) );
  memset( &refVec, 0, sizeof(refVec) );
  memset( &xVec, 0, sizeof(xVec) );
  memset( &yVec, 0, sizeof(yVec) );
    {
        char * fname;
        const char* dir;
        dir=context->isFull ? FULL_VECTOR_DIR :BRIEF_VECTOR_DIR;
        fname=(char*)malloc(strlen(context_fft->fInName)+strlen(dir)+2);
        sprintf(fname,"%s/%s",dir,context_fft->fInName);
        fIn = fopen(fname, "rb");
        free(fname);
        fname=(char*)malloc(strlen(context_fft->fRefName)+strlen(dir)+2);
        sprintf(fname,"%s/%s",dir,context_fft->fRefName);
        fRef = fopen(fname, "rb");
        free(fname);
    }

  /* Allocate vectors for in/out/reference data. */
  if ( !vecAlloc( &inVec , lenIn , TE_ALIGN_YES, FMT_FRACT16|fmtX, 0 ) ||
       !vecAlloc( &outVec, lenOut, TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) ||
       !vecAlloc( &refVec, N     , TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) )
  {
    printf( "fileProc_stnd(): failed to allocate in/out/ref vectors, N=%d\n", N );
  }
  /* Open input data file. */
  else if (fIn==NULL )
  {
    printf( "fileProc_stnd(): failed to open %s for reading\n", context_fft->fInName );
  }
  /* Open reference data file. */
  else if (fRef == NULL )
  {
    printf( "fileProc_stnd(): failed to open %s for reading\n", context_fft->fRefName );
  }
  else res = 1;

  /*
   * Process the input file frame-by-frame.
   */

  if ( res )
  {
    int baseFmt = ( context->desc->fmt & FMT_DTYPE_MASK );
    int isAligned = context->desc->isAligned;

    int16_t   * pin  = (int16_t   *)vecGetElem( &inVec , 0 );
    float64_t * pout = (float64_t *)vecGetElem( &outVec, 0 );
    float64_t * pref = (float64_t *)vecGetElem( &refVec, 0 );

    float32_t sinadAvg, sinadMin = HUGE_VALF;
    float64_t errSum = 0, refSum = 0;

    int efbMin = 32;

    memset( &xVec, 0, sizeof(xVec) );
    memset( &yVec, 0, sizeof(yVec) );

    context_fft->frameCnt = 0;

    frameProcFxn = ( (tTestEngFrameProcFxn_fft_stnd*)context->target.fut );

    /* Read 16-bit complex samples from the input file. */
    while ((n = fread(pin, inVec.szElem, lenFread, fIn)) > 0)
    {
      /* Augment the (last) incomplete frame with zeros. */
        memset((uint8_t*)pin + n*inVec.szElem, 0, (lenFread - n)*inVec.szElem);
      /* Zero the output frame. */
      memset( pout, 0, lenOut*outVec.szElem );

      /* For inverse FFT test, unpack the Nyquist frequency bin. */
      if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
      {

        pin[2*(lenIn-1)+0] = pin[1];
        pin[2*(lenIn-1)+1] = pin[1] = 0;
      }

      /* Allocate in/out buffers for the target FFT routine. */
      if ( ( !xVec.szBulk && !vecAlloc( &xVec, lenIn , isAligned, baseFmt | fmtX, 0 ) ) ||
           ( !yVec.szBulk && !vecAlloc( &yVec, lenOut, isAligned, baseFmt | fmtY, 0 ) ) )
      {
        printf( "fileProc_stnd(): failed to allocate xVec/yVec, "
                "frameCnt=%d, N=%d\n", context_fft->frameCnt, N );
        res = 0; break;
      }
	  context_fft->dataSize = (size_t)lenIn+(size_t)lenOut;

      /* Use a proprietary frame processing function to perform the target FFT. */
      if ( !( res = frameProcFxn( context, pin, pout, &xVec, &yVec ) ) ) break;

      /* When in unaligned mode, in/out vectors should be periodically reallocated to try various
       * address offsets. */
      if ( !( ++context_fft->frameCnt % XY_VEC_REALLOC_PERIOD ) && !isAligned )
      {
        if ( !( res = vecsFree( &xVec, &yVec, 0 ) ) )
        {
          printf( "fileProc_stnd(): vecsFree() failed for xVec/yVec, "
                  "frameCnt=%d, N=%d\n", context_fft->frameCnt, N );
          break;
        }

        memset( &xVec, 0, sizeof(xVec) );
        memset( &yVec, 0, sizeof(yVec) );
      }

#ifdef PSEUDOFLOAT
	  static float64_t output_fl32[2*8192];
	  static float64_t output_fl64[2*8192];

	  float32_fft(pin, output_fl32, N);
	  float64_fft(pin, output_fl64, N);	  
#endif

      /* Read double precision reference samples. Reference frame always contains
       * exactly N samples (real or complex)! */
      if ( (int)fread( pref, refVec.szElem, N, fRef ) < N )
      {
        printf( "fileProc_stnd(): failed to read reference data, "
                "frameCnt=%d, N=%d\n", context_fft->frameCnt, N );
        res = 0; break;
      }

      /* Estimate the SINAD ratio for the current frame, and update error/reference
       * power sums for the whole file SINAD. */
      {
        float64_t refMax = 0, errMax = 0;
        float64_t p, q;
        int m, M;

        M = ( ( fmtY & FMT_CPLX ) ? 2*lenOut : lenOut );

        for ( p=q=0, m=0; m<M; m++ )
        {
          float64_t err = pout[m] - pref[m];

          /* |signal|^2 */
          p += pref[m]*pref[m];
          /* |noise+distortion|^2 */
          q += err*err;

          refMax = MAX( refMax, fabs( pref[m] ) );
          errMax = MAX( errMax, fabs( err     ) );
        }

        refSum += p;
        errSum += q;

        if ( p>0 )
        {
          sinadMin = MIN( sinadMin, (float32_t)(p/q) );
          efbMin   = MAX( 0, MIN( efbMin, L_sub_ll( (int)logb(refMax), (int)logb(errMax) ) ) );
        }
      }
    }

    /*
     * Finalize the min SINAD estimation and print a summary.
     */

    if ( res )
    {
      sinadMin = 10.f*log10f( sinadMin );

      /* Set the test result for test verification. */
      *vecGetElem_fl32( &context->dataSet.Z, 0 ) = sinadMin;

      if ( context->isVerbose )
      {
        sinadAvg = ( refSum > 0 ? 10.f*log10f( (float32_t)(refSum/errSum) ) : -HUGE_VALF );

        if ( ( baseFmt & FMT_DTYPE_MASK ) == FMT_FRACT16 || ( baseFmt & FMT_DTYPE_MASK ) == FMT_FRACT32 )
        {
          printf( "Error-Free Bits %2d  SINAD min %5.1f dB  avg %5.1f dB  ",
                  efbMin, sinadMin, sinadAvg );
        }
        else
        {
          printf( "SINAD min %5.1f dB  avg %5.1f dB  ", sinadMin, sinadAvg );
        }
      }
    }
  }

  /*
   * Close files and free vectors.
   */

  if ( fIn  ) fclose( fIn  );
  if ( fRef ) fclose( fRef );

  if ( inVec .szBulk > 0 ) vecFree( &inVec  );
  if ( outVec.szBulk > 0 ) vecFree( &outVec );
  if ( refVec.szBulk > 0 ) vecFree( &refVec );
  if ( xVec  .szBulk > 0 ) vecFree( &xVec   );
  if ( yVec  .szBulk > 0 ) vecFree( &yVec   );

} /* fileProc_stnd() */

#ifdef PSEUDOFLOAT
#define PI 3.14159265
void float32_fft(int16_t *in, float64_t *out, int N)
{
	int k, n;
	float64_t angle;
	float32_t sum;

	for (k = 0; k < N; k++)
	{
		out[2 * k] = 0;
		out[2 * k + 1] = 0;
		sum = 0;
		for (n = 0; n < N; n++)
		{
			angle = 2.0 * PI*k*n / (float32_t)N;
			sum += in[n] * cos(angle);
		}
		out[2 * k]     = ldexp(sum, -15);

		sum = 0;
		for (n = 0; n < N; n++) 
		{
			angle = 2.0 * PI*k*n / (float32_t)N;
			sum -= in[n] * sin(angle);
		}
		out[2 * k + 1] = ldexp(sum, -15);;
	}
}

void float64_fft(int16_t *in, float64_t *out, int N)
{
	int k, n;
	float64_t multVal, angle;

	for (k = 0; k < (N / 2) + 1; k++)
	{
		out[2 * k] = 0;
		out[2 * k + 1] = 0;
		for (n = 0; n < N; n++)
		{
			angle = (float64_t)2.0 * PI*k*n / (float64_t)N;
			multVal = ((float64_t) in[n]) * (float64_t)cos(angle);
			out[2 * k] += multVal;
		}
		out[2 * k] = (double) (ldexp(out[2 * k], -15));

		for (n = 0; n < N; n++)
		{
			angle = (float64_t)2.0 * PI*k*n / (float64_t)N;
			multVal = ((float64_t)in[n]) * (float64_t)sin(angle);
			out[2 * k + 1] -= multVal;
		}
		out[2 * k + 1] = (double)(ldexp(out[2 * k + 1], -15));
	}
}

void float32_ifft(int16_t* in, float64_t* out, int N)
{            
    int k, n;
    float64_t angle;
    float32_t sum;
            
    for (k = 0; k < N; k++)
    {       
        out[k] = 0;
        sum = 0;
        for (n = 0; n < N; n++)
        {   
            angle = 2.0 * PI * k * n / (float32_t)N;
            sum = sum + in[n] * cos(angle) - in[n + 1] * sin(angle);
            // X[n]  = X[n] + X_Real[k] * cos((2 * M_PI * k * n) / N) - X_Imag[k] * sin((2 * M_PI * k * n) / N);            
        }   
        out[k] = ldexp(sum, -15);
    }   
}           

void float64_ifft(int16_t* in, float64_t* out, int N)
{
    int k, n;
    float64_t sum, angle;

    for (k = 0; k < N ; k++)
    {
        out[k] = 0;
        sum = 0;
        for (n = 0; n < N; n++)
        {
            angle = (float64_t)2.0 * PI * k * n / (float64_t)N;
            sum = sum + ((float64_t)in[n]) * cos(angle) - ((float64_t)in[n]) * sin(angle);
        }
        out[k] = (double)(ldexp(sum, -15));
    }
}



#endif