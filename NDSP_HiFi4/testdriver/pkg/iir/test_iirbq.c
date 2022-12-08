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
 * Test procedures for IIR
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(iir)
/* Test engine extension for IIR filters. */
#include "testeng_iir.h"
#include "testeng_iir_old.h"

#if 0 //HiFi3/3z API
static size_t bqriir24x24_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR24X24_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir24x24_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR24X24_DF2_SCRATCH_SIZE( N, M ) ; }
static const tIirOldDescr api_bqriir24x24_df1={(tIirOldFxnAlloc*)bqriir24x24_df1_alloc,(tIirOldFxnInit*)bqriir24x24_df1_init,bqriir24x24_df1_getScratch,(tIirOldFxnProcess*)bqriir24x24_df1};
static const tIirOldDescr api_bqriir24x24_df2={(tIirOldFxnAlloc*)bqriir24x24_df2_alloc,(tIirOldFxnInit*)bqriir24x24_df2_init,bqriir24x24_df2_getScratch,(tIirOldFxnProcess*)bqriir24x24_df2};
static const tTestEngDesc descr_bqriir24x24_df1 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1,                NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_old, te_destroy_iir_old, &te_loadFxn_iir_old, &te_processFxn_iir_old };
static const tTestEngDesc descr_bqriir24x24_df2 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF2,                NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_old, te_destroy_iir_old, &te_loadFxn_iir_old, &te_processFxn_iir_old };
#endif
static size_t bqriir16x16_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir16x16_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR16X16_DF2_SCRATCH_SIZE( N, M ) ; }

static size_t bqriir32x16_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir32x16_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X16_DF2_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir32x32_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir32x32_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X32_DF2_SCRATCH_SIZE( N, M ) ; }
static size_t stereo_bqriir16x16_df1_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t stereo_bqriir32x16_df1_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR32X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t stereo_bqriir32x32_df1_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) ; }

static size_t bqriir16x16_df1_nd_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR16X16_DF1_ND_SCRATCH_SIZE(N, M); }
static size_t bqriir16x16_df2_nd_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR16X16_DF2_ND_SCRATCH_SIZE(N, M); }

static size_t bqriir32x16_df1_nd_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X16_DF1_ND_SCRATCH_SIZE(N, M); }
static size_t bqriir32x16_df2_nd_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X16_DF2_ND_SCRATCH_SIZE(N, M); }
static size_t bqriir32x32_df1_nd_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X32_DF1_ND_SCRATCH_SIZE(N, M); }
static size_t bqriir32x32_df2_nd_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X32_DF2_ND_SCRATCH_SIZE(N, M); }

static size_t stereo_bqriir16x16_df1_nd_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR16X16_DF1_ND_SCRATCH_SIZE(N, M); }
static size_t stereo_bqriir32x16_df1_nd_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR32X16_DF1_ND_SCRATCH_SIZE(N, M); }
static size_t stereo_bqriir32x32_df1_nd_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR32X32_DF1_ND_SCRATCH_SIZE(N, M); }

static const tIirDescr api_bqriirf_df1 = { (tIirFxnAlloc*)bqriirf_df1_alloc, (tIirFxnInit*)bqriirf_df1_init, NULL, (tIirFxnProcess*)bqriirf_df1, (tIirFxnDelay*)bqriirf_df1_groupDelay };
static const tIirDescr api_bqriirf_df2 = { (tIirFxnAlloc*)bqriirf_df2_alloc, (tIirFxnInit*)bqriirf_df2_init, NULL, (tIirFxnProcess*)bqriirf_df2, (tIirFxnDelay*)bqriirf_df2_groupDelay };
static const tIirDescr api_bqriirf_df2t = { (tIirFxnAlloc*)bqriirf_df2t_alloc, (tIirFxnInit*)bqriirf_df2t_init, NULL, (tIirFxnProcess*)bqriirf_df2t, (tIirFxnDelay*)bqriirf_df2t_groupDelay };
static const tIirDescr api_bqciirf_df1 = { (tIirFxnAlloc*)bqciirf_df1_alloc, (tIirFxnInit*)bqciirf_df1_init, NULL, (tIirFxnProcess*)bqciirf_df1, (tIirFxnDelay*)bqciirf_df1_groupDelay };
static const tIirDescr api_bqriir16x16_df1 = { (tIirFxnAlloc*)bqriir16x16_df1_alloc, (tIirFxnInit*)bqriir16x16_df1_init, bqriir16x16_df1_getScratch, (tIirFxnProcess*)bqriir16x16_df1, (tIirFxnDelay*)bqriir16x16_df1_groupDelay };
static const tIirDescr api_bqriir16x16_df2 = { (tIirFxnAlloc*)bqriir16x16_df2_alloc, (tIirFxnInit*)bqriir16x16_df2_init, bqriir16x16_df2_getScratch, (tIirFxnProcess*)bqriir16x16_df2, (tIirFxnDelay*)bqriir16x16_df2_groupDelay };
static const tIirDescr api_bqriir32x16_df1 = { (tIirFxnAlloc*)bqriir32x16_df1_alloc, (tIirFxnInit*)bqriir32x16_df1_init, bqriir32x16_df1_getScratch, (tIirFxnProcess*)bqriir32x16_df1, (tIirFxnDelay*)bqriir32x16_df1_groupDelay };
static const tIirDescr api_bqriir32x16_df2 = { (tIirFxnAlloc*)bqriir32x16_df2_alloc, (tIirFxnInit*)bqriir32x16_df2_init, bqriir32x16_df2_getScratch, (tIirFxnProcess*)bqriir32x16_df2, (tIirFxnDelay*)bqriir32x16_df2_groupDelay };
static const tIirDescr api_bqriir32x32_df1 = { (tIirFxnAlloc*)bqriir32x32_df1_alloc, (tIirFxnInit*)bqriir32x32_df1_init, bqriir32x32_df1_getScratch, (tIirFxnProcess*)bqriir32x32_df1, (tIirFxnDelay*)bqriir32x32_df1_groupDelay };
static const tIirDescr api_bqriir32x32_df2 = { (tIirFxnAlloc*)bqriir32x32_df2_alloc, (tIirFxnInit*)bqriir32x32_df2_init, bqriir32x32_df2_getScratch, (tIirFxnProcess*)bqriir32x32_df2, (tIirFxnDelay*)bqriir32x32_df2_groupDelay };
static const tIirStereoDescr api_stereo_bqriirf_df1 = { (tIirStereoFxnAlloc*)stereo_bqriirf_df1_alloc, (tIirStereoFxnInit*)stereo_bqriirf_df1_init, NULL, (tIirStereoFxnProcess*)stereo_bqriirf_df1, (tIirStereoFxnDelay*)stereo_bqriirf_df1_groupDelay };
static const tIirStereoDescr api_stereo_bqriir16x16_df1 = { (tIirStereoFxnAlloc*)stereo_bqriir16x16_df1_alloc, (tIirStereoFxnInit*)stereo_bqriir16x16_df1_init, stereo_bqriir16x16_df1_getScratch, (tIirStereoFxnProcess*)stereo_bqriir16x16_df1, (tIirStereoFxnDelay*)stereo_bqriir16x16_df1_groupDelay };
static const tIirStereoDescr api_stereo_bqriir32x16_df1 = { (tIirStereoFxnAlloc*)stereo_bqriir32x16_df1_alloc, (tIirStereoFxnInit*)stereo_bqriir32x16_df1_init, stereo_bqriir32x16_df1_getScratch, (tIirStereoFxnProcess*)stereo_bqriir32x16_df1, (tIirStereoFxnDelay*)stereo_bqriir32x16_df1_groupDelay };
static const tIirStereoDescr api_stereo_bqriir32x32_df1 = { (tIirStereoFxnAlloc*)stereo_bqriir32x32_df1_alloc, (tIirStereoFxnInit*)stereo_bqriir32x32_df1_init, stereo_bqriir32x32_df1_getScratch, (tIirStereoFxnProcess*)stereo_bqriir32x32_df1, (tIirStereoFxnDelay*)stereo_bqriir32x32_df1_groupDelay };


static const tIirDescr api_bqriirf_df1_nd = { (tIirFxnAlloc*)bqriirf_df1_nd_alloc, (tIirFxnInit*)bqriirf_df1_nd_init, NULL, (tIirFxnProcess*)bqriirf_df1_nd, (tIirFxnDelay*)bqriirf_df1_nd_groupDelay };
static const tIirDescr api_bqriirf_df2_nd = { (tIirFxnAlloc*)bqriirf_df2_nd_alloc, (tIirFxnInit*)bqriirf_df2_nd_init, NULL, (tIirFxnProcess*)bqriirf_df2_nd, (tIirFxnDelay*)bqriirf_df2_nd_groupDelay };
static const tIirDescr api_bqriirf_df2t_nd = { (tIirFxnAlloc*)bqriirf_df2t_nd_alloc, (tIirFxnInit*)bqriirf_df2t_nd_init, NULL, (tIirFxnProcess*)bqriirf_df2t_nd, (tIirFxnDelay*)bqriirf_df2t_nd_groupDelay };
static const tIirDescr api_bqciirf_df1_nd = { (tIirFxnAlloc*)bqciirf_df1_nd_alloc, (tIirFxnInit*)bqciirf_df1_nd_init, NULL, (tIirFxnProcess*)bqciirf_df1_nd, (tIirFxnDelay*)bqciirf_df1_nd_groupDelay };
static const tIirDescr api_bqriir16x16_df1_nd = { (tIirFxnAlloc*)bqriir16x16_df1_nd_alloc, (tIirFxnInit*)bqriir16x16_df1_nd_init, bqriir16x16_df1_nd_getScratch, (tIirFxnProcess*)bqriir16x16_df1_nd, (tIirFxnDelay*)bqriir16x16_df1_nd_groupDelay };
static const tIirDescr api_bqriir16x16_df2_nd = { (tIirFxnAlloc*)bqriir16x16_df2_nd_alloc, (tIirFxnInit*)bqriir16x16_df2_nd_init, bqriir16x16_df2_nd_getScratch, (tIirFxnProcess*)bqriir16x16_df2_nd, (tIirFxnDelay*)bqriir16x16_df2_nd_groupDelay };
static const tIirDescr api_bqriir32x16_df1_nd = { (tIirFxnAlloc*)bqriir32x16_df1_nd_alloc, (tIirFxnInit*)bqriir32x16_df1_nd_init, bqriir32x16_df1_nd_getScratch, (tIirFxnProcess*)bqriir32x16_df1_nd, (tIirFxnDelay*)bqriir32x16_df1_nd_groupDelay };
static const tIirDescr api_bqriir32x16_df2_nd = { (tIirFxnAlloc*)bqriir32x16_df2_nd_alloc, (tIirFxnInit*)bqriir32x16_df2_nd_init, bqriir32x16_df2_nd_getScratch, (tIirFxnProcess*)bqriir32x16_df2_nd, (tIirFxnDelay*)bqriir32x16_df2_nd_groupDelay };
static const tIirDescr api_bqriir32x32_df1_nd = { (tIirFxnAlloc*)bqriir32x32_df1_nd_alloc, (tIirFxnInit*)bqriir32x32_df1_nd_init, bqriir32x32_df1_nd_getScratch, (tIirFxnProcess*)bqriir32x32_df1_nd, (tIirFxnDelay*)bqriir32x32_df1_nd_groupDelay };
static const tIirDescr api_bqriir32x32_df2_nd = { (tIirFxnAlloc*)bqriir32x32_df2_nd_alloc, (tIirFxnInit*)bqriir32x32_df2_nd_init, bqriir32x32_df2_nd_getScratch, (tIirFxnProcess*)bqriir32x32_df2_nd, (tIirFxnDelay*)bqriir32x32_df2_nd_groupDelay };
static const tIirStereoDescr api_stereo_bqriirf_df1_nd = { (tIirStereoFxnAlloc*)stereo_bqriirf_df1_nd_alloc, (tIirStereoFxnInit*)stereo_bqriirf_df1_nd_init, NULL, (tIirStereoFxnProcess*)stereo_bqriirf_df1_nd, (tIirStereoFxnDelay*)stereo_bqriirf_df1_nd_groupDelay };
static const tIirStereoDescr api_stereo_bqriir16x16_df1_nd = { (tIirStereoFxnAlloc*)stereo_bqriir16x16_df1_nd_alloc, (tIirStereoFxnInit*)stereo_bqriir16x16_df1_nd_init, stereo_bqriir16x16_df1_nd_getScratch, (tIirStereoFxnProcess*)stereo_bqriir16x16_df1_nd, (tIirStereoFxnDelay*)stereo_bqriir16x16_df1_nd_groupDelay };
static const tIirStereoDescr api_stereo_bqriir32x16_df1_nd = { (tIirStereoFxnAlloc*)stereo_bqriir32x16_df1_nd_alloc, (tIirStereoFxnInit*)stereo_bqriir32x16_df1_nd_init, stereo_bqriir32x16_df1_nd_getScratch, (tIirStereoFxnProcess*)stereo_bqriir32x16_df1_nd, (tIirStereoFxnDelay*)stereo_bqriir32x16_df1_nd_groupDelay };
static const tIirStereoDescr api_stereo_bqriir32x32_df1_nd = { (tIirStereoFxnAlloc*)stereo_bqriir32x32_df1_nd_alloc, (tIirStereoFxnInit*)stereo_bqriir32x32_df1_nd_init, stereo_bqriir32x32_df1_nd_getScratch, (tIirStereoFxnProcess*)stereo_bqriir32x32_df1_nd, (tIirStereoFxnDelay*)stereo_bqriir32x32_df1_nd_groupDelay };

static const tTestEngDesc descr_bqriirf_df1 = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2 = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2t = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2T | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqciirf_df1 = { FMT_CPLX | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir16x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir16x16_df2 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF2 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x16_df2 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF2, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x32_df1 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x32_df2 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF2, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_stereo_bqriirf_df1 = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir16x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir32x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir32x32_df1 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };

static const tTestEngDesc descr_bqriirf_df1_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2t_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2T | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqciirf_df1_nd = { FMT_CPLX | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir16x16_df1_nd = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir16x16_df2_nd = { FMT_REAL | FMT_FRACT16, TE_IIR_DF2 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x16_df1_nd = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x16_df2_nd = { FMT_REAL | FMT_FRACT16, TE_IIR_DF2, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x32_df1_nd = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x32_df2_nd = { FMT_REAL | FMT_FRACT32, TE_IIR_DF2, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_stereo_bqriirf_df1_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir16x16_df1_nd = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir32x16_df1_nd = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir32x32_df1_nd = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };



typedef struct
{
  int                 phaseNum;
  const tTestEngDesc *pIirDescr;
  tTestEngTarget      fxns;
  int                 runAlways;   /* 1 - brief & full, 0 - full only */
  const char*         seqFile;
}
tTbl;

static const tTbl tests[] =
{
    //tests with delay imitation 
    { 1, &descr_bqriir16x16_df1, (tTestEngTarget)&api_bqriir16x16_df1, 1, "bqriir16x16_df1_bpf1.seq" },
    { 1, &descr_bqriir16x16_df1, (tTestEngTarget)&api_bqriir16x16_df1, 1, "bqriir16x16_df1_lpf1.seq" },    
    { 1, &descr_bqriir16x16_df1, (tTestEngTarget)&api_bqriir16x16_df1, 0, "bqriir16x16_df1_bpf2.seq" },
    { 1, &descr_bqriir16x16_df1, (tTestEngTarget)&api_bqriir16x16_df1, 0, "bqriir16x16_df1_bsf1.seq" },
    { 1, &descr_bqriir16x16_df1, (tTestEngTarget)&api_bqriir16x16_df1, 0, "bqriir16x16_df1_hpf1.seq" },

    { 1, &descr_bqriir16x16_df2, (tTestEngTarget)&api_bqriir16x16_df2, 1, "bqriir16x16_df2_lpf1.seq" },
    { 1, &descr_bqriir16x16_df2, (tTestEngTarget)&api_bqriir16x16_df2, 1, "bqriir16x16_df2_bpf1.seq" },
    { 1, &descr_bqriir16x16_df2, (tTestEngTarget)&api_bqriir16x16_df2, 0, "bqriir16x16_df2_bpf2.seq" },
  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 0, "bqriir16x16_df2_bsf1.seq" },
  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 0, "bqriir16x16_df2_hpf1.seq" },

  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 1, "bqriir32x16_df1_lpf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 1, "bqriir32x16_df1_bpf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "bqriir32x16_df1_bpf2.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "bqriir32x16_df1_bsf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "bqriir32x16_df1_hpf1.seq" },

  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 1, "bqriir32x16_df2_lpf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 1, "bqriir32x16_df2_bpf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "bqriir32x16_df2_bpf2.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "bqriir32x16_df2_bsf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "bqriir32x16_df2_hpf1.seq" },
#if 0 //HiFi3/3z API
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 1, "bqriir24x24_df1_lpf1.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 1, "bqriir24x24_df1_bpf1.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 0, "bqriir24x24_df1_bpf2.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 0, "bqriir24x24_df1_bsf1.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 0, "bqriir24x24_df1_hpf1.seq" },

  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 1, "bqriir24x24_df2_lpf1.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 1, "bqriir24x24_df2_bpf1.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 0, "bqriir24x24_df2_bpf2.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 0, "bqriir24x24_df2_bsf1.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 0, "bqriir24x24_df2_hpf1.seq" },
#endif
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 1, "bqriir32x32_df1_lpf1.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 1, "bqriir32x32_df1_bpf1.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "bqriir32x32_df1_bpf2.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "bqriir32x32_df1_bsf1.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "bqriir32x32_df1_hpf1.seq" },

  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 1, "bqriir32x32_df2_lpf1.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 1, "bqriir32x32_df2_bpf1.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "bqriir32x32_df2_bpf2.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "bqriir32x32_df2_bsf1.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "bqriir32x32_df2_hpf1.seq" },

  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 1, "stereo_bqriir16x16_df1_lpf1.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 1, "stereo_bqriir16x16_df1_bpf1.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 0, "stereo_bqriir16x16_df1_bpf2.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 0, "stereo_bqriir16x16_df1_bsf1.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 0, "stereo_bqriir16x16_df1_hpf1.seq" },

  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 1, "stereo_bqriir32x16_df1_lpf1.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 1, "stereo_bqriir32x16_df1_bpf1.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 0, "stereo_bqriir32x16_df1_bpf2.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 0, "stereo_bqriir32x16_df1_bsf1.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 0, "stereo_bqriir32x16_df1_hpf1.seq" },

  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 1, "stereo_bqriir32x32_df1_lpf1.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 1, "stereo_bqriir32x32_df1_bpf1.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 0, "stereo_bqriir32x32_df1_bpf2.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 0, "stereo_bqriir32x32_df1_bsf1.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 0, "stereo_bqriir32x32_df1_hpf1.seq" },

  /*
   * Stage 2
   */
  { 2, &descr_bqriirf_df1    , (tTestEngTarget)&api_bqriirf_df1    , 1, "bqriirf_df1_lpf1.seq" },
  { 2, &descr_bqriirf_df1    , (tTestEngTarget)&api_bqriirf_df1    , 1, "bqriirf_df1_bpf1.seq" },
  { 2, &descr_bqriirf_df1    , (tTestEngTarget)&api_bqriirf_df1    , 0, "bqriirf_df1_bpf2.seq" },
  { 2, &descr_bqriirf_df1    , (tTestEngTarget)&api_bqriirf_df1    , 0, "bqriirf_df1_bpf3.seq" },
  { 2, &descr_bqriirf_df1    , (tTestEngTarget)&api_bqriirf_df1    , 0, "bqriirf_df1_bsf1.seq" },
  { 2, &descr_bqriirf_df1    , (tTestEngTarget)&api_bqriirf_df1    , 0, "bqriirf_df1_hpf1.seq" },

  { 2, &descr_bqriirf_df2    , (tTestEngTarget)&api_bqriirf_df2    , 1, "bqriirf_df2_lpf1.seq" },
  { 2, &descr_bqriirf_df2    , (tTestEngTarget)&api_bqriirf_df2    , 1, "bqriirf_df2_bpf1.seq" },
  { 2, &descr_bqriirf_df2    , (tTestEngTarget)&api_bqriirf_df2    , 0, "bqriirf_df2_bpf2.seq" },
  { 2, &descr_bqriirf_df2    , (tTestEngTarget)&api_bqriirf_df2    , 0, "bqriirf_df2_bpf3.seq" },
  { 2, &descr_bqriirf_df2    , (tTestEngTarget)&api_bqriirf_df2    , 0, "bqriirf_df2_bsf1.seq" },
  { 2, &descr_bqriirf_df2    , (tTestEngTarget)&api_bqriirf_df2    , 0, "bqriirf_df2_hpf1.seq" },

  { 2, &descr_bqriirf_df2t   , (tTestEngTarget)&api_bqriirf_df2t   , 1, "bqriirf_df2t_lpf1.seq" },
  { 2, &descr_bqriirf_df2t   , (tTestEngTarget)&api_bqriirf_df2t   , 1, "bqriirf_df2t_bpf1.seq" },
  { 2, &descr_bqriirf_df2t   , (tTestEngTarget)&api_bqriirf_df2t   , 0, "bqriirf_df2t_bpf2.seq" },
  { 2, &descr_bqriirf_df2t   , (tTestEngTarget)&api_bqriirf_df2t   , 0, "bqriirf_df2t_bpf3.seq" },
  { 2, &descr_bqriirf_df2t   , (tTestEngTarget)&api_bqriirf_df2t   , 0, "bqriirf_df2t_bsf1.seq" },
  { 2, &descr_bqriirf_df2t   , (tTestEngTarget)&api_bqriirf_df2t   , 0, "bqriirf_df2t_hpf1.seq" },

  { 2, &descr_bqciirf_df1    , (tTestEngTarget)&api_bqciirf_df1    , 1, "bqciirf_df1_lpf1.seq" },
  { 2, &descr_bqciirf_df1    , (tTestEngTarget)&api_bqciirf_df1    , 1, "bqciirf_df1_bpf1.seq" },
  { 2, &descr_bqciirf_df1    , (tTestEngTarget)&api_bqciirf_df1    , 0, "bqciirf_df1_bpf2.seq" },
  { 2, &descr_bqciirf_df1    , (tTestEngTarget)&api_bqciirf_df1    , 0, "bqciirf_df1_bsf1.seq" },
  { 2, &descr_bqciirf_df1    , (tTestEngTarget)&api_bqciirf_df1    , 0, "bqciirf_df1_hpf1.seq" },

  { 2, &descr_stereo_bqriirf_df1    , (tTestEngTarget)&api_stereo_bqriirf_df1    , 1, "stereo_bqriirf_df1_lpf1.seq" },
  { 2, &descr_stereo_bqriirf_df1    , (tTestEngTarget)&api_stereo_bqriirf_df1    , 1, "stereo_bqriirf_df1_bpf1.seq" },
  { 2, &descr_stereo_bqriirf_df1    , (tTestEngTarget)&api_stereo_bqriirf_df1    , 0, "stereo_bqriirf_df1_bpf2.seq" },
  { 2, &descr_stereo_bqriirf_df1    , (tTestEngTarget)&api_stereo_bqriirf_df1    , 0, "stereo_bqriirf_df1_bpf3.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 0, "stereo_bqriirf_df1_bsf1.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 0, "stereo_bqriirf_df1_hpf1.seq" },

};

static const tTbl tests_nd[] =
{
    //no delay functions
    { 1, &descr_bqriir16x16_df1_nd, (tTestEngTarget)&api_bqriir16x16_df1_nd, 1, "bqriir16x16_df1_nd_bpf1.seq" },
    { 1, &descr_bqriir16x16_df1_nd, (tTestEngTarget)&api_bqriir16x16_df1_nd, 1, "bqriir16x16_df1_nd_lpf1.seq" },
    { 1, &descr_bqriir16x16_df1_nd, (tTestEngTarget)&api_bqriir16x16_df1_nd, 0, "bqriir16x16_df1_nd_bpf2.seq" },
    { 1, &descr_bqriir16x16_df1_nd, (tTestEngTarget)&api_bqriir16x16_df1_nd, 0, "bqriir16x16_df1_nd_bsf1.seq" },
    { 1, &descr_bqriir16x16_df1_nd, (tTestEngTarget)&api_bqriir16x16_df1_nd, 0, "bqriir16x16_df1_nd_hpf1.seq" },

    { 1, &descr_bqriir16x16_df2_nd, (tTestEngTarget)&api_bqriir16x16_df2_nd, 1, "bqriir16x16_df2_nd_lpf1.seq" },
    { 1, &descr_bqriir16x16_df2_nd, (tTestEngTarget)&api_bqriir16x16_df2_nd, 1, "bqriir16x16_df2_nd_bpf1.seq" },
    { 1, &descr_bqriir16x16_df2_nd, (tTestEngTarget)&api_bqriir16x16_df2_nd, 0, "bqriir16x16_df2_nd_bpf2.seq" },
    { 1, &descr_bqriir16x16_df2_nd, (tTestEngTarget)&api_bqriir16x16_df2_nd, 0, "bqriir16x16_df2_nd_bsf1.seq" },
    { 1, &descr_bqriir16x16_df2_nd, (tTestEngTarget)&api_bqriir16x16_df2_nd, 0, "bqriir16x16_df2_nd_hpf1.seq" },

    { 1, &descr_bqriir32x16_df1_nd, (tTestEngTarget)&api_bqriir32x16_df1_nd, 1, "bqriir32x16_df1_nd_lpf1.seq" },
    { 1, &descr_bqriir32x16_df1_nd, (tTestEngTarget)&api_bqriir32x16_df1_nd, 1, "bqriir32x16_df1_nd_bpf1.seq" },
    { 1, &descr_bqriir32x16_df1_nd, (tTestEngTarget)&api_bqriir32x16_df1_nd, 0, "bqriir32x16_df1_nd_bpf2.seq" },
    { 1, &descr_bqriir32x16_df1_nd, (tTestEngTarget)&api_bqriir32x16_df1_nd, 0, "bqriir32x16_df1_nd_bsf1.seq" },
    { 1, &descr_bqriir32x16_df1_nd, (tTestEngTarget)&api_bqriir32x16_df1_nd, 0, "bqriir32x16_df1_nd_hpf1.seq" },

    { 1, &descr_bqriir32x16_df2_nd, (tTestEngTarget)&api_bqriir32x16_df2_nd, 1, "bqriir32x16_df2_nd_lpf1.seq" },
    { 1, &descr_bqriir32x16_df2_nd, (tTestEngTarget)&api_bqriir32x16_df2_nd, 1, "bqriir32x16_df2_nd_bpf1.seq" },
    { 1, &descr_bqriir32x16_df2_nd, (tTestEngTarget)&api_bqriir32x16_df2_nd, 0, "bqriir32x16_df2_nd_bpf2.seq" },
    { 1, &descr_bqriir32x16_df2_nd, (tTestEngTarget)&api_bqriir32x16_df2_nd, 0, "bqriir32x16_df2_nd_bsf1.seq" },
    { 1, &descr_bqriir32x16_df2_nd, (tTestEngTarget)&api_bqriir32x16_df2_nd, 0, "bqriir32x16_df2_nd_hpf1.seq" },

    { 1, &descr_bqriir32x32_df1_nd, (tTestEngTarget)&api_bqriir32x32_df1_nd, 1, "bqriir32x32_df1_nd_lpf1.seq" },
    { 1, &descr_bqriir32x32_df1_nd, (tTestEngTarget)&api_bqriir32x32_df1_nd, 1, "bqriir32x32_df1_nd_bpf1.seq" },
    { 1, &descr_bqriir32x32_df1_nd, (tTestEngTarget)&api_bqriir32x32_df1_nd, 0, "bqriir32x32_df1_nd_bpf2.seq" },
    { 1, &descr_bqriir32x32_df1_nd, (tTestEngTarget)&api_bqriir32x32_df1_nd, 0, "bqriir32x32_df1_nd_bsf1.seq" },
    { 1, &descr_bqriir32x32_df1_nd, (tTestEngTarget)&api_bqriir32x32_df1_nd, 0, "bqriir32x32_df1_nd_hpf1.seq" },

    { 1, &descr_bqriir32x32_df2_nd, (tTestEngTarget)&api_bqriir32x32_df2_nd, 1, "bqriir32x32_df2_nd_lpf1.seq" },
    { 1, &descr_bqriir32x32_df2_nd, (tTestEngTarget)&api_bqriir32x32_df2_nd, 1, "bqriir32x32_df2_nd_bpf1.seq" },
    { 1, &descr_bqriir32x32_df2_nd, (tTestEngTarget)&api_bqriir32x32_df2_nd, 0, "bqriir32x32_df2_nd_bpf2.seq" },
    { 1, &descr_bqriir32x32_df2_nd, (tTestEngTarget)&api_bqriir32x32_df2_nd, 0, "bqriir32x32_df2_nd_bsf1.seq" },
    { 1, &descr_bqriir32x32_df2_nd, (tTestEngTarget)&api_bqriir32x32_df2_nd, 0, "bqriir32x32_df2_nd_hpf1.seq" },

    { 1, &descr_stereo_bqriir16x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir16x16_df1_nd, 1, "stereo_bqriir16x16_df1_nd_lpf1.seq" },
    { 1, &descr_stereo_bqriir16x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir16x16_df1_nd, 1, "stereo_bqriir16x16_df1_nd_bpf1.seq" },
    { 1, &descr_stereo_bqriir16x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir16x16_df1_nd, 0, "stereo_bqriir16x16_df1_nd_bpf2.seq" },
    { 1, &descr_stereo_bqriir16x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir16x16_df1_nd, 0, "stereo_bqriir16x16_df1_nd_bsf1.seq" },
    { 1, &descr_stereo_bqriir16x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir16x16_df1_nd, 0, "stereo_bqriir16x16_df1_nd_hpf1.seq" },

    { 1, &descr_stereo_bqriir32x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x16_df1_nd, 1, "stereo_bqriir32x16_df1_nd_lpf1.seq" },
    { 1, &descr_stereo_bqriir32x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x16_df1_nd, 1, "stereo_bqriir32x16_df1_nd_bpf1.seq" },
    { 1, &descr_stereo_bqriir32x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x16_df1_nd, 0, "stereo_bqriir32x16_df1_nd_bpf2.seq" },
    { 1, &descr_stereo_bqriir32x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x16_df1_nd, 0, "stereo_bqriir32x16_df1_nd_bsf1.seq" },
    { 1, &descr_stereo_bqriir32x16_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x16_df1_nd, 0, "stereo_bqriir32x16_df1_nd_hpf1.seq" },

    { 1, &descr_stereo_bqriir32x32_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x32_df1_nd, 1, "stereo_bqriir32x32_df1_nd_lpf1.seq" },
    { 1, &descr_stereo_bqriir32x32_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x32_df1_nd, 1, "stereo_bqriir32x32_df1_nd_bpf1.seq" },
    { 1, &descr_stereo_bqriir32x32_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x32_df1_nd, 0, "stereo_bqriir32x32_df1_nd_bpf2.seq" },
    { 1, &descr_stereo_bqriir32x32_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x32_df1_nd, 0, "stereo_bqriir32x32_df1_nd_bsf1.seq" },
    { 1, &descr_stereo_bqriir32x32_df1_nd, (tTestEngTarget)&api_stereo_bqriir32x32_df1_nd, 0, "stereo_bqriir32x32_df1_nd_hpf1.seq" },

    /*
    * Stage 2
    */
    //no delay functions
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 1, "bqriirf_df1_nd_lpf1.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 1, "bqriirf_df1_nd_bpf1.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "bqriirf_df1_nd_bpf2.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "bqriirf_df1_nd_bpf3.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "bqriirf_df1_nd_bsf1.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "bqriirf_df1_nd_hpf1.seq" },

    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 1, "bqriirf_df2_nd_lpf1.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 1, "bqriirf_df2_nd_bpf1.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "bqriirf_df2_nd_bpf2.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "bqriirf_df2_nd_bpf3.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "bqriirf_df2_nd_bsf1.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "bqriirf_df2_nd_hpf1.seq" },

    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 1, "bqriirf_df2t_nd_lpf1.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 1, "bqriirf_df2t_nd_bpf1.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "bqriirf_df2t_nd_bpf2.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "bqriirf_df2t_nd_bpf3.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "bqriirf_df2t_nd_bsf1.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "bqriirf_df2t_nd_hpf1.seq" },

    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 1, "bqciirf_df1_nd_lpf1.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 1, "bqciirf_df1_nd_bpf1.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 0, "bqciirf_df1_nd_bpf2.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 0, "bqciirf_df1_nd_bsf1.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 0, "bqciirf_df1_nd_hpf1.seq" },

    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 1, "stereo_bqriirf_df1_nd_lpf1.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 1, "stereo_bqriirf_df1_nd_bpf1.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "stereo_bqriirf_df1_nd_bpf2.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "stereo_bqriirf_df1_nd_bpf3.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "stereo_bqriirf_df1_nd_bsf1.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "stereo_bqriirf_df1_nd_hpf1.seq" },

};



/* Perform all tests for IIR API functions. */
int main_iirbq( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
    int n;
    int res = 1;
    for (n=0; n<(int)(sizeof(tests)/sizeof(tests[0])); n++)
    {
        if ( ( phaseNum == 0 || phaseNum == tests[n].phaseNum ) && ( isFull || tests[n].runAlways ) )
        {
            res &= (0!=TestEngRun(tests[n].fxns, tests[n].pIirDescr, tests[n].seqFile, isFull, isVerbose, breakOnError,0));
            if (res == 0 && breakOnError) break;
        }
    }
    return (res);
} /* main_iirbq() */

int main_iirbq_nd(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int n;
    int res = 1;
    for (n = 0; n<(int)(sizeof(tests_nd) / sizeof(tests_nd[0])); n++)
    {
        if ((phaseNum == 0 || phaseNum == tests_nd[n].phaseNum) && (isFull || tests_nd[n].runAlways))
        {
            res &= (0 != TestEngRun(tests_nd[n].fxns, tests_nd[n].pIirDescr, tests_nd[n].seqFile, isFull, isVerbose, breakOnError, 0));
            if (res == 0 && breakOnError) break;
        }
    }
    return (res);
} /* main_iirbq_nd() */

