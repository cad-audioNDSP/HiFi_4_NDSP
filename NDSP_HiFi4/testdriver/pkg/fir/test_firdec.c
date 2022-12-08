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

#include <string.h>
#include <stdlib.h>
#if 0// for HiFi3/3z
static tFirOldDescr api_firdec24x24     ={{(tFirOldFxnAlloc*)firdec24x24_alloc,     (tFirOldFxnInit*)firdec24x24_init,     (tFirOldFxnProcess*)firdec24x24_process     }};
static const tTestEngDesc descr_firdec24x24     = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
#endif

/* phase 1*/
static tFirOldDescr api_firdec16x16     ={{(tFirOldFxnAlloc*)firdec16x16_alloc,     (tFirOldFxnInit*)firdec16x16_init,     (tFirOldFxnProcess*)firdec16x16_process     }};
static tFirOldDescr api_firdec32x16     ={{(tFirOldFxnAlloc*)firdec32x16_alloc,     (tFirOldFxnInit*)firdec32x16_init,     (tFirOldFxnProcess*)firdec32x16_process     }};
static tFirOldDescr api_firdec32x32     ={{(tFirOldFxnAlloc*)firdec32x32_alloc,     (tFirOldFxnInit*)firdec32x32_init,     (tFirOldFxnProcess*)firdec32x32_process     }};
static tFirOldDescr api_firdec32x32ep   ={{(tFirOldFxnAlloc*)firdec32x32ep_alloc,   (tFirOldFxnInit*)firdec32x32ep_init,   (tFirOldFxnProcess*)firdec32x32ep_process   }};

static const tTestEngDesc descr_firdec16x16     = { FMT_REAL|FMT_FRACT16, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firdec32x16     = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR|TE_FIR_FILTER_32X16, NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firdec32x32     = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firdec32x32ep   = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
/* phase 2*/
static tFirOldDescr api_firdecf       ={{(tFirOldFxnAlloc*)firdecf_alloc,       (tFirOldFxnInit*)firdecf_init,       (tFirOldFxnProcess*)firdecf_process       }};
static const tTestEngDesc descr_firdecf         = { FMT_REAL|FMT_FLOAT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR             ,NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

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
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec16x16_dn2x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec16x16_dn3x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec16x16_dn11x.seq"},
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec16x16_dn6x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,0,"firdec16x16_dn4x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,0,"firdec16x16_dn5x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,0,"firdec16x16_dn23x.seq"},

  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec32x16_dn2x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec32x16_dn3x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec32x16_dn11x.seq"},
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec32x16_dn6x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,0,"firdec32x16_dn4x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,0,"firdec32x16_dn5x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,0,"firdec32x16_dn23x.seq"},
#if 0// for HiFi3/3z
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn2x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn3x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn11x.seq"},
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn6x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,0,"firdec24x24_dn4x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,0,"firdec24x24_dn5x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,0,"firdec24x24_dn23x.seq"},
#endif
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec32x32_dn2x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec32x32_dn3x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec32x32_dn11x.seq"},
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec32x32_dn6x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,0,"firdec32x32_dn4x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,0,"firdec32x32_dn5x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,0,"firdec32x32_dn23x.seq"},

  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec32x32ep_dn2x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec32x32ep_dn3x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec32x32ep_dn11x.seq"},
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec32x32ep_dn6x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,0,"firdec32x32ep_dn4x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,0,"firdec32x32ep_dn5x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,0,"firdec32x32ep_dn23x.seq"},

  /*
   * Stage 2
   */

  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdecf_dn2x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdecf_dn3x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdecf_dn11x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdecf_dn6x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,0,"firdecf_dn4x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,0,"firdecf_dn5x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,0,"firdecf_dn23x.seq"},

};

/* Perform all tests for FIR API functions. */
int main_firdec( int phaseNum, int isFull, int isVerbose, int breakOnError )
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
} /* main_firdec() */

