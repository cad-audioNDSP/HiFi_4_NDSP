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
#include "testeng_fir_old.h"

#if 0// for HiFi3/3z
static tFirOldDescr api_firinterp24x24  ={{(tFirOldFxnAlloc*)firinterp24x24_alloc,  (tFirOldFxnInit*)firinterp24x24_init,  (tFirOldFxnProcess*)firinterp24x24_process  }};
static const tTestEngDesc descr_firinterp24x24  = { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
#endif

/* phase 1*/
static tFirOldDescr api_firinterp16x16  ={{(tFirOldFxnAlloc*)firinterp16x16_alloc,  (tFirOldFxnInit*)firinterp16x16_init,  (tFirOldFxnProcess*)firinterp16x16_process  }};
static tFirOldDescr api_firinterp32x16  ={{(tFirOldFxnAlloc*)firinterp32x16_alloc,  (tFirOldFxnInit*)firinterp32x16_init,  (tFirOldFxnProcess*)firinterp32x16_process  }};
static tFirOldDescr api_firinterp32x32  ={{(tFirOldFxnAlloc*)firinterp32x32_alloc,  (tFirOldFxnInit*)firinterp32x32_init,  (tFirOldFxnProcess*)firinterp32x32_process  }};
static tFirOldDescr api_firinterp32x32ep={{(tFirOldFxnAlloc*)firinterp32x32ep_alloc,(tFirOldFxnInit*)firinterp32x32ep_init,(tFirOldFxnProcess*)firinterp32x32ep_process}};

static const tTestEngDesc descr_firinterp16x16  = { FMT_REAL | FMT_FRACT16, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firinterp32x16  = { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR|TE_FIR_FILTER_32X16, NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firinterp32x32  = { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firinterp32x32ep= { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
/* phase 2*/
static tFirOldDescr api_firinterpf    ={{(tFirOldFxnAlloc*)firinterpf_alloc,    (tFirOldFxnInit*)firinterpf_init,    (tFirOldFxnProcess*)firinterpf_process    }};
static const tTestEngDesc descr_firinterpf      = { FMT_REAL | FMT_FLOAT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR             ,NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

typedef struct
{
  int                 phaseNum;
  const tTestEngDesc *pFirDescr;
  tTestEngTarget      fxns;
  int                 runAlways;   /* 1 - brief & full, 0 - full only */
  const char*         seqFile;
}
tTbl;

static const tTbl tests[] =
{
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,1,"firinterp16x16_up2x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,1,"firinterp16x16_up3x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,1,"firinterp16x16_up6x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,0,"firinterp16x16_up4x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,0,"firinterp16x16_up5x.seq" },

  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,1,"firinterp32x16_up2x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,1,"firinterp32x16_up3x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,1,"firinterp32x16_up6x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,0,"firinterp32x16_up4x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,0,"firinterp32x16_up5x.seq" },
#if 0// for HiFi3/3z
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,1,"firinterp24x24_up2x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,1,"firinterp24x24_up3x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,1,"firinterp24x24_up6x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,0,"firinterp24x24_up4x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,0,"firinterp24x24_up5x.seq" },
#endif
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,1,"firinterp32x32_up2x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,1,"firinterp32x32_up3x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,1,"firinterp32x32_up6x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,0,"firinterp32x32_up4x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,0,"firinterp32x32_up5x.seq" },

  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,1,"firinterp32x32ep_up2x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,1,"firinterp32x32ep_up3x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,1,"firinterp32x32ep_up6x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,0,"firinterp32x32ep_up4x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,0,"firinterp32x32ep_up5x.seq" },

  /*
   * Stage 2
   */
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,1,"firinterpf_up2x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,1,"firinterpf_up3x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,1,"firinterpf_up6x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,0,"firinterpf_up4x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,0,"firinterpf_up5x.seq"},
};

/* Perform all tests for FIR API functions. */
int main_firint( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
    int n;
    int res = 1;
    for (n=0; n<(int)(sizeof(tests)/sizeof(tests[0])); n++)
    {
        if ( ( phaseNum == 0 || phaseNum == tests[n].phaseNum ) && ( isFull || tests[n].runAlways ) )
        {
            res &= (0!=TestEngRun(tests[n].fxns, tests[n].pFirDescr, tests[n].seqFile, isFull, isVerbose, breakOnError,0));
            if (res == 0 && breakOnError) break;
        }
    }

    return (res);
} /* main_firint() */

