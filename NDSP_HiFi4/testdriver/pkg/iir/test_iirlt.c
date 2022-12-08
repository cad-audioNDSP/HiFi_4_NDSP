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
/* Test engine extension for latiice IIR filters. */
#include "testeng_iir_lat.h"

#if 0 //HiFi3/3z API
static const tIirLatDescr api_latr24x24={(tIirLatFxnAlloc*)latr24x24_alloc,(tIirLatFxnInit*)latr24x24_init,(tIirLatFxnProcess*)latr24x24_process};
static const tTestEngDesc descr_latr24x24     = { FMT_REAL | FMT_FRACT32, 0, NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_lat, te_destroy_iir_lat, &te_loadFxn_iir_lat, &te_processFxn_iir_lat };
#endif
static const tIirLatDescr api_latr16x16={(tIirLatFxnAlloc*)latr16x16_alloc,(tIirLatFxnInit*)latr16x16_init,(tIirLatFxnProcess*)latr16x16_process};
static const tIirLatDescr api_latr32x16={(tIirLatFxnAlloc*)latr32x16_alloc,(tIirLatFxnInit*)latr32x16_init,(tIirLatFxnProcess*)latr32x16_process};
static const tIirLatDescr api_latr32x32={(tIirLatFxnAlloc*)latr32x32_alloc,(tIirLatFxnInit*)latr32x32_init,(tIirLatFxnProcess*)latr32x32_process};
static const tIirLatDescr api_latrf    ={(tIirLatFxnAlloc*)latrf_alloc,    (tIirLatFxnInit*)latrf_init,    (tIirLatFxnProcess*)latrf_process    };
static const tTestEngDesc descr_latr16x16     = { FMT_REAL | FMT_FRACT16, TE_IIR_16X16, NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_lat, te_destroy_iir_lat, &te_loadFxn_iir_lat, &te_processFxn_iir_lat };
static const tTestEngDesc descr_latr32x16     = { FMT_REAL | FMT_FRACT16, 0, NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_lat, te_destroy_iir_lat, &te_loadFxn_iir_lat, &te_processFxn_iir_lat };
static const tTestEngDesc descr_latr32x32     = { FMT_REAL | FMT_FRACT32, 0, NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_lat, te_destroy_iir_lat, &te_loadFxn_iir_lat, &te_processFxn_iir_lat };
static const tTestEngDesc descr_latrf         = { FMT_REAL | FMT_FLOAT32, 0, NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_lat, te_destroy_iir_lat, &te_loadFxn_iir_lat, &te_processFxn_iir_lat };


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
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 1, "latr16x16_lpf1.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 0, "latr16x16_lpf2.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 0, "latr16x16_lpf3.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 0, "latr16x16_lpf4.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 1, "latr16x16_lpf5.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 0, "latr16x16_lpf6.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 0, "latr16x16_lpf7.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 0, "latr16x16_lpf8.seq"       },
  { 1, &descr_latr16x16          , (tTestEngTarget)&api_latr16x16          , 1, "latr16x16_lpf9.seq"       },
#if 0 //HiFi3/3z API
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 1, "latr24x24_lpf1.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 0, "latr24x24_lpf2.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 0, "latr24x24_lpf3.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 0, "latr24x24_lpf4.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 1, "latr24x24_lpf5.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 0, "latr24x24_lpf6.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 0, "latr24x24_lpf7.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 0, "latr24x24_lpf8.seq"       },
  { 1, &descr_latr24x24          , (tTestEngTarget)&api_latr24x24          , 1, "latr24x24_lpf9.seq"       },
#endif                                                                                                          
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 1, "latr32x16_lpf1.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 0, "latr32x16_lpf2.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 0, "latr32x16_lpf3.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 0, "latr32x16_lpf4.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 1, "latr32x16_lpf5.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 0, "latr32x16_lpf6.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 0, "latr32x16_lpf7.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 0, "latr32x16_lpf8.seq"       },
  { 1, &descr_latr32x16          , (tTestEngTarget)&api_latr32x16          , 1, "latr32x16_lpf9.seq"       },

  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 1, "latr32x32_lpf1.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 0, "latr32x32_lpf2.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 0, "latr32x32_lpf3.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 0, "latr32x32_lpf4.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 1, "latr32x32_lpf5.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 0, "latr32x32_lpf6.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 0, "latr32x32_lpf7.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 0, "latr32x32_lpf8.seq"       },
  { 1, &descr_latr32x32          , (tTestEngTarget)&api_latr32x32          , 1, "latr32x32_lpf9.seq"       },

  /*
   * Stage 2
   */

  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 1, "latrf_lpf1.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "latrf_lpf2.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "latrf_lpf3.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "latrf_lpf4.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 1, "latrf_lpf5.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "latrf_lpf6.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "latrf_lpf7.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "latrf_lpf8.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 1, "latrf_lpf9.seq"           },
};


/* Perform all tests for IIR API functions. */
int main_iirlt( int phaseNum, int isFull, int isVerbose, int breakOnError )
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
} /* main_iirlt() */
