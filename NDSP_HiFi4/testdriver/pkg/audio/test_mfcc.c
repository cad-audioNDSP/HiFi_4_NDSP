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
 * Test procesdures for MFCC features extraction APIs.
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(audio)
/* Test engine API. */
#include "testeng.h"
/* Test engine add-on for log mel filterbank and MFCC features extractor tests. */
#include "testeng_logmel.h"
#include "testeng_mfcc.h"

static const te_logmel_api_t logmel32x32_api = {
    (te_logmel_alloc_fxn_t          *)&logmel32x32_alloc,
    (te_logmel_init_fxn_t           *)&logmel32x32_init,
    (te_logmel_process_fxn_t        *)&logmel32x32_process,
    (te_logmel_getScratchSize_fxn_t *)&logmel32x32_getScratchSize,
};

static const te_logmel_api_t logmelf_api = {
    (te_logmel_alloc_fxn_t          *)&logmelf_alloc,
    (te_logmel_init_fxn_t           *)&logmelf_init,
    (te_logmel_process_fxn_t        *)&logmelf_process,
    (te_logmel_getScratchSize_fxn_t *)&logmelf_getScratchSize,
};

static const te_mfcc_api_t mfcc32x32_api = {
    (te_mfcc_getDefaultParams_fxn_t *)&mfcc_getDefaultParams,
    (te_mfcc_alloc_fxn_t            *)&mfcc32x32_alloc,
    (te_mfcc_init_fxn_t             *)&mfcc32x32_init,
    (te_mfcc_process_fxn_t          *)&mfcc32x32_process,
    (te_mfcc_getScratchSize_fxn_t   *)&mfcc32x32_getScratchSize
};

static const te_mfcc_api_t mfccf_api = {
    (te_mfcc_getDefaultParams_fxn_t *)&mfcc_getDefaultParams,
    (te_mfcc_alloc_fxn_t            *)&mfccf_alloc,
    (te_mfcc_init_fxn_t             *)&mfccf_init,
    (te_mfcc_process_fxn_t          *)&mfccf_process,
    (te_mfcc_getScratchSize_fxn_t   *)&mfccf_getScratchSize
};

#define LOGMEL_TEST_DESC(fmt, align) { (fmt), 0, NULL, TE_DIM_NUM_1, (align), &te_createFxn_logmel, &te_destroyFxn_logmel, \
                                                                              &te_loadFxn_logmel, &te_processFxn_logmel }
#define MFCC_TEST_DESC(fmt, align)   { (fmt), 0, NULL, TE_DIM_NUM_2, (align), &te_createFxn_mfcc, &te_destroyFxn_mfcc, \
                                                                              &te_loadFxn_mfcc, &te_processFxn_mfcc }

static const tTestEngDesc logmel32x32_desc = LOGMEL_TEST_DESC(FMT_REAL|FMT_FRACT32, TE_ALIGN_YES);
static const tTestEngDesc logmelf_desc     = LOGMEL_TEST_DESC(FMT_REAL|FMT_FLOAT32, TE_ALIGN_YES);
static const tTestEngDesc mfcc32x32_desc   = MFCC_TEST_DESC(FMT_REAL|FMT_FRACT32, TE_ALIGN_YES);
static const tTestEngDesc mfccf_desc       = MFCC_TEST_DESC(FMT_REAL|FMT_FLOAT32, TE_ALIGN_YES);

#define DO_TEST( api, desc, seqFile )                                                          \
    if ( res || !breakOnError ) res &= ( 0 != TestEngRun((tTestEngTarget)&api, &desc, seqFile, \
                                                         isFull, isVerbose, breakOnError, 0) )

/* Perform all tests for MFCC features extraction APIs. */
int main_mfcc( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
    int res = 1;
    if (phaseNum==0 || phaseNum==1) {
        DO_TEST(logmel32x32_api, logmel32x32_desc, "logmel32x32.seq");
        DO_TEST(  mfcc32x32_api,   mfcc32x32_desc,   "mfcc32x32.seq");
    }
    if (phaseNum==0 || phaseNum==2) {
        DO_TEST(logmelf_api, logmelf_desc, "logmelf.seq");
        DO_TEST(  mfccf_api,   mfccf_desc,   "mfccf.seq");
    }
    return (res);
} /* main_mfcc() */
