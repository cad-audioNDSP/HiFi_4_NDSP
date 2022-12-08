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
* Test procedures for complex FFT functions.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

#if !defined(COMPILER_MSVC)
/* C99 compatible type */
#include <complex.h>
#endif

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Suppress Visual C warnings on +/-HUGE_VALF macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

#define sz_fr16c    sizeof(complex_fract16)
#define sz_fl64c    sizeof(complex_double)

/* Period, in frames, between reallocation of in/out buffers consumed by the
* target FFT routine. */
#define XY_VEC_REALLOC_PERIOD     16

/* Test data file processing for a standalone test. Applies the target FFT function
* to test data, compares the FFT output with reference data and estimates the SINAD. */
static void fileProc_stnd(tTestEngContext * context);


/* Initializer for a test description structure. */
#define TEST_DESC_STND( fmt, extraParam, align,fwdfft,invfft )  { { (fmt),(extraParam), NULL, TE_ARGNUM_1,(align), \
    &te_create_fft, &te_destroy_fft, \
    &te_load_fft, &fileProc_stnd }, (tFrwTransFxn)fwdfft, (tInvTransFxn)invfft}

/* FFT attributes  */
#define ATTR_FAST               (TE_FFT_CPLX|TE_FFT_OPT_REUSE_INBUF)

typedef const struct tagTestDef
{
    int                phaseNum;
    const char *       seqFilename;
    tTestEngTarget     target;
    tTestEngDesc_fft   desc;

} TestDef_t;

TestDef_t testTbl_cfft[] =
{
#if 1
    { 1, "fft16_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft32_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft64_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft128_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft256_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft512_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft1024_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft2048_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft4096_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },

    { 1, "fft16_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft32_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft64_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft128_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft256_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft512_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft1024_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft2048_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "fft4096_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },

    { 1, "ifft16_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft32_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft64_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft128_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft256_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft512_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft1024_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft2048_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft4096_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },

    { 1, "ifft16_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft32_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft64_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft128_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft256_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft512_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft1024_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft2048_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "ifft4096_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
#endif
#if 0 //HiFi3/3z API
    { 1, "fft_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "fft_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "fft_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "fft_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "ifft_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "ifft_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "ifft_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "ifft_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
#endif

#if 1
    { 1, "fft16_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft32_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft64_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft128_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft256_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft512_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft1024_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft2048_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft4096_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft16_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft32_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft64_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft128_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft256_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft512_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft1024_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft2048_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "fft4096_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },

    { 1, "ifft16_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft32_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft64_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft128_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft256_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft512_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft1024_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft2048_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft4096_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft16_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft32_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft64_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft128_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft256_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft512_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft1024_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft2048_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "ifft4096_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
#endif

#if 1
    { 1, "fft16_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft32_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft64_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft128_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft256_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft512_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft1024_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft2048_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft4096_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },

    { 1, "fft16_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft32_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft64_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft128_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft256_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft512_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft1024_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft2048_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft4096_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },

    { 1, "ifft16_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft32_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft64_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft128_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft256_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft512_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft1024_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft2048_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft4096_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },

    { 1, "ifft16_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft32_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft64_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft128_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft256_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft512_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft1024_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft2048_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft4096_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
#endif


    { 0 } /* End of table */
}; // testTbl_cfft

TestDef_t testTbl_cfftie[] =
{
#if 1
    { 1, "fft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },
    { 1, "fft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },
    { 1, "fft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },
    { 1, "fft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },

    { 1, "ifft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
    { 1, "ifft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
    { 1, "ifft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
    { 1, "ifft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
#endif
#if 0 //HiFi3/3z API
    { 1, "fft_cplx24x24_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx24x24_ie, NULL) },
    { 1, "ifft_cplx24x24_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx24x24_ie) },
#endif
#if 1
    { 1, "fft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "fft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "fft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },

    { 1, "fft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "fft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "fft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },

    { 1, "ifft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "ifft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "ifft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },

    { 1, "ifft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "ifft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "ifft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
#endif
#if 1
    { 1, "fft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "fft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "fft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "fft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },

    { 1, "fft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "fft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "fft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "fft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },

    { 1, "ifft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "ifft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "ifft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "ifft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },

    { 1, "ifft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "ifft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "ifft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "ifft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
#endif
#if 1
    { 2, "fft8_cplxf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft16_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft32_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft64_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft128_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft256_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft512_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft1024_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft2048_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },
    { 2, "fft4096_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, &fft_cplxf_ie, NULL) },

    { 2, "ifft8_cplxf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft16_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft32_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft64_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft128_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft256_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft512_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft1024_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft2048_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
    { 2, "ifft4096_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_cplxf_ie) },
#endif
#if 1
    { 1, "stereo_fft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },
    { 1, "stereo_fft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },
    { 1, "stereo_fft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },
    { 1, "stereo_fft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },

    { 1, "stereo_ifft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
    { 1, "stereo_ifft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
    { 1, "stereo_ifft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
    { 1, "stereo_ifft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
#endif
#if 1
    { 1, "stereo_fft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "stereo_fft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "stereo_fft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },

    { 1, "stereo_fft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "stereo_fft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "stereo_fft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },

    { 1, "stereo_ifft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "stereo_ifft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "stereo_ifft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },

    { 1, "stereo_ifft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "stereo_ifft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "stereo_ifft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
#endif
#if 1
    { 1, "stereo_fft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "stereo_fft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "stereo_fft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "stereo_fft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },

    { 1, "stereo_fft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "stereo_fft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "stereo_fft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "stereo_fft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },

    { 1, "stereo_ifft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "stereo_ifft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "stereo_ifft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "stereo_ifft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },

    { 1, "stereo_ifft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "stereo_ifft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "stereo_ifft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "stereo_ifft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
#endif
#if 1
    { 2, "stereo_fft8_cplxf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft16_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft32_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft64_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft128_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft256_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft512_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft1024_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft2048_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },
    { 2, "stereo_fft4096_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplxf_ie, NULL) },

    { 2, "stereo_ifft8_cplxf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft16_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft32_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft64_cplxf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft128_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft256_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft512_cplxf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft1024_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft2048_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
    { 2, "stereo_ifft4096_cplxf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplxf_ie) },
#endif
    { 0 } /* End of table */
}; // testTbl_cfftie

TestDef_t testTbl_cnfft[] =
{
#if 1
    // Mixed radix forward 16x16
    { 1, "fft160_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft192_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft240_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft320_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft384_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft480_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft160_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft192_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft240_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft320_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft384_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    { 1, "fft480_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx16x16, NULL) },
    // Mixed radix inverse 16x16
    { 1, "ifft160_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft192_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft240_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft320_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft384_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft480_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft160_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft192_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft240_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft320_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft384_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
    { 1, "ifft480_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx16x16) },
#endif
#if 1
    // Mixed radix forward 32x16
    { 1, "fft160_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft192_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft240_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft320_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft384_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft480_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft160_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft192_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft240_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft320_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft384_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    { 1, "fft480_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x16, NULL) },
    // Mixed radix inverse 32x16
    { 1, "ifft160_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft192_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft240_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft320_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft384_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft480_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft160_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft192_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft240_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft320_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft384_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
    { 1, "ifft480_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_32X16 | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x16) },
#endif
#if 1
    // Mixed radix forward 32x32
    { 1, "fft12_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft24_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft36_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft48_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft60_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft72_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft80_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft96_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft100_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft108_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft120_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft144_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft160_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft180_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft192_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft200_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft216_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft240_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft288_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft300_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft320_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft324_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft360_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft384_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft400_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft432_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft480_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft540_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft576_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft600_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft768_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft960_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft12_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft24_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft36_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft48_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft60_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft72_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft80_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft96_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft100_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft108_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft120_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft144_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft160_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft180_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft192_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft200_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft216_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft240_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft288_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft300_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft320_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft324_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft360_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft384_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft400_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft432_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft480_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft540_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft576_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft600_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft768_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "fft960_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    // Mixed radix inverse 32x32
    { 1, "ifft12_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft24_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft36_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft48_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft60_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft72_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft80_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft96_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft100_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft108_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft120_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft144_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft160_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft180_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft192_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft200_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft216_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft240_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft288_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft300_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft320_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft324_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft360_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft384_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft400_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft432_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft480_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft540_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft576_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft600_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft768_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft960_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft12_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft24_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft36_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft48_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft60_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft72_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft80_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft96_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft100_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft108_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft120_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft144_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft160_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft180_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft192_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft200_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft216_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft240_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft288_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft300_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft320_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft324_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft360_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft384_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft400_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft432_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft480_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft540_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft576_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft600_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft768_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "ifft960_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
#endif

    { 0 } /* End of table */
}; // testTbl_cfftie

/* Perform all tests for fft_cplx, ifft_cplx API functions. */
int main_cfft(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_cfft[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_cfft[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_cfft[ix].target;

            if (!IS_PRESENT(testTbl_cfft[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_cfft[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_cfft[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_cfft[ix].desc.desc,
                testTbl_cfft[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);

} /* main_cfft() */


/* Perform all tests for fft_cplx*_ie, ifft_cplx*_ie API functions. */
int main_cfftie(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_cfftie[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_cfftie[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_cfftie[ix].target;

            if (!IS_PRESENT(testTbl_cfftie[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_cfftie[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_cfftie[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_cfftie[ix].desc.desc,
                testTbl_cfftie[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);

} /* main_cfftie() */

/* Perform all tests for mixed radix fft_cplx32x32, ifft_cplx32x32 API functions. */
int main_cnfft(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_cnfft[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_cnfft[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_cnfft[ix].target;

            if (!IS_PRESENT(testTbl_cnfft[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_cnfft[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_cnfft[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_cnfft[ix].desc.desc,
                testTbl_cnfft[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);

} /* main_cnfft() */


/* Test data file processing for a standalone test. Applies the target FFT function
* to test data, compares the FFT output with reference data and estimates the SINAD. */
void fileProc_stnd(tTestEngContext * context)
{
    tTestEngContext_fft           * context_fft;
    tTestEngFrameProcFxn_fft_stnd * frameProcFxn;

    tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
    tVec xVec, yVec;            /* In/out vectors for the target FFT routine. */
    FILE *fIn = 0, *fRef = 0;
    int res = 0;
    int n, N;

    NASSERT(context && context->desc);
    NASSERT(context->target.fut && context->target.handle);

    context_fft = (tTestEngContext_fft*)context->target.handle;

    /* FFT size */
    N = context->args.N_DIM_;
    if (context->desc->extraParam & TE_FFT_STEREO) N*=2;

    /* Preset an invalid SINAD for the test result. */
    *vecGetElem_fl32(&context->dataSet.Z, 0) = -HUGE_VALF;

    if (context->isVerbose)
    {
        if (context->desc->extraParam & TE_FFT_OPT_SCALE_METH)
        {
            printf("scale_mtd %d ", context_fft->scale_method);
        }

        printf("%-40s  ", context_fft->fRefName);
    }

    memset(&inVec, 0, sizeof(inVec));
    memset(&outVec, 0, sizeof(outVec));
    memset(&refVec, 0, sizeof(refVec));
    memset(&xVec, 0, sizeof(xVec));
    memset(&yVec, 0, sizeof(yVec));
    {
        char * fname;
        const char *dir;        
        dir=context->isFull ? FULL_VECTOR_DIR :BRIEF_VECTOR_DIR;
        fname=(char *)malloc(strlen(context_fft->fInName)+strlen(dir)+2);
        sprintf(fname,"%s/%s",dir,context_fft->fInName);
        fIn = fopen(fname, "rb");
        free(fname);
        fname=(char *)malloc(strlen(context_fft->fRefName)+strlen(dir)+2);
        sprintf(fname,"%s/%s",dir,context_fft->fRefName);
        fRef = fopen(fname, "rb");
        free(fname);
    }

    /* Allocate vectors for in/out/reference data. */
    if (!vecAlloc(&inVec, N, TE_ALIGN_YES, FMT_FRACT16 | FMT_CPLX, 0) ||
        !vecAlloc(&outVec, N, TE_ALIGN_YES, FMT_FLOAT64 | FMT_CPLX, 0) ||
        !vecAlloc(&refVec, N, TE_ALIGN_YES, FMT_FLOAT64 | FMT_CPLX, 0))
    {
        printf("fileProc_stnd(): failed to allocate in/out/ref vectors, N=%d\n", N);
    }
    /* Open input data file. */
    else if (fIn==NULL)
    {
        printf("fileProc_stnd(): failed to open %s for reading\n", context_fft->fInName);
    }
    /* Open reference data file. */
    else if (fRef==NULL)
    {
        printf("fileProc_stnd(): failed to open %s for reading\n", context_fft->fRefName);
    }
    else res = 1;

    /*
    * Process the input file frame-by-frame.
    */

    if (res)
    {
        int fmt = context->desc->fmt | FMT_CPLX;
        int isAligned = context->desc->isAligned;

        complex_fract16 * pin = vecGetElem_fr16c(&inVec, 0);
        complex_double  * pout = vecGetElem_fl64c(&outVec, 0);
        complex_double  * pref = vecGetElem_fl64c(&refVec, 0);

        float32_t sinadAvg, sinadMin = HUGE_VALF;
        float64_t errSum = 0, refSum = 0;

        int efbMin = 32;

        memset(&xVec, 0, sizeof(xVec));
        memset(&yVec, 0, sizeof(yVec));

        context_fft->frameCnt = 0;

        frameProcFxn = ((tTestEngFrameProcFxn_fft_stnd*)context->target.fut);

        /* Read N 16-bit complex samples from the input file. */
        while ((n = fread(pin, sz_fr16c, N, fIn)) > 0)
        {
            /* Augment the (last) incomplete frame with zeros. */
            memset((uint8_t*)pin + n*sz_fr16c, 0, (N - n)*sz_fr16c);
            /* Zero the output frame. */
            memset(pout, 0, N*sz_fl64c);

            /* Allocate in/out buffers for the target FFT routine. */
            if (!xVec.szBulk && (2 != vecsAlloc(isAligned, fmt, &xVec, N, &yVec, N, 0)))
            {
                printf("fileProc_stnd(): failed to allocate xVec/yVec, "
                    "frameCnt=%d, fmt=0x%x, N=%d\n", context_fft->frameCnt, fmt, N);
                res = 0; break;
            }
            context_fft->dataSize = (size_t)xVec.szElem*xVec.nElem + (size_t)yVec.szElem*yVec.nElem;

            /* Use a proprietary frame processing function to perform the target FFT. */
            if (!(res = frameProcFxn(context, (int16_t*)pin, (float64_t*)pout, &xVec, &yVec))) break;

            /* When in unaligned mode, in/out vectors should be periodically reallocated to try various
            * address offsets. */
            if (!(++context_fft->frameCnt % XY_VEC_REALLOC_PERIOD) && !isAligned)
            {
                if (!(res = vecsFree(&xVec, &yVec, 0)))
                {
                    printf("fileProc_stnd(): vecsFree() failed for xVec/yVec, "
                        "frameCnt=%d, fmt=0x%x, N=%d\n", context_fft->frameCnt, fmt, N);
                    break;
                }

                memset(&xVec, 0, sizeof(xVec));
                memset(&yVec, 0, sizeof(yVec));
            }

            /* Read N double precision complex reference samples. */
            if ((int)fread(pref, sz_fl64c, N, fRef) < N)
            {
                printf("fileProc_stnd(): failed to read reference data, "
                    "frameCnt=%d, fmt=0x%x, N=%d\n", context_fft->frameCnt, fmt, N);
                res = 0; break;
            }

            /* Estimate the SINAD ratio for the current frame, and update error/reference
            * power sums for the whole file SINAD. */
            {
                float64_t refMax = 0, errMax = 0;
                float64_t p, q;

                for (p = q = 0, n = 0; n<N; n++)
                {
                    double err_re, err_im;
                    err_re = creal(pout[n]) - creal(pref[n]);
                    err_im = cimag(pout[n]) - cimag(pref[n]);

                    /* |signal|^2 */
                    p += creal(pref[n])*creal(pref[n]) + cimag(pref[n])*cimag(pref[n]);
                    /* |noise+distortion|^2 */
                    q += err_re*err_re + err_im*err_im;

                    refMax = MAX(refMax, MAX(fabs(creal(pref[n])), fabs(cimag(pref[n]))));
                    errMax = MAX(errMax, MAX(fabs(err_re), fabs(err_im)));
                }

                refSum += p;
                errSum += q;

                if (p>0)
                {
                    sinadMin = MIN(sinadMin, (float32_t)(p / q));
                    efbMin = MAX(0, MIN(efbMin, L_sub_ll((int)logb(refMax), (int)logb(errMax))));
                }
            }
        }

        /*
        * Finalize the min SINAD estimation and print a summary.
        */

        if (res)
        {
            sinadMin = 10.f*log10f(sinadMin);

            /* Set the test result for test verification. */
            *vecGetElem_fl32(&context->dataSet.Z, 0) = sinadMin;

            if (context->isVerbose)
            {
                sinadAvg = (refSum > 0 ? 10.f*log10f((float32_t)(refSum / errSum)) : -HUGE_VALF);

                if ((fmt & FMT_DTYPE_MASK) == FMT_FRACT16 || (fmt & FMT_DTYPE_MASK) == FMT_FRACT32)
                {
                    printf("Error-Free Bits %2d  SINAD min %5.1f dB  avg %5.1f dB  ",
                        efbMin, sinadMin, sinadAvg);
                }
                else
                {
                    printf("SINAD min %5.1f dB  avg %5.1f dB  ", sinadMin, sinadAvg);
                }
            }
        }
    }

    /*
    * Close files and free vectors.
    */

    if (fIn) fclose(fIn);
    if (fRef) fclose(fRef);

    if (inVec.szBulk > 0) vecFree(&inVec);
    if (outVec.szBulk > 0) vecFree(&outVec);
    if (refVec.szBulk > 0) vecFree(&refVec);
    if (xVec.szBulk > 0) vecFree(&xVec); //!!! xVec not initialized !!!
    if (yVec.szBulk > 0) vecFree(&yVec);

} /* fileProc_stnd() */
