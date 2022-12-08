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
 * Test procedures for matrix functions
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matop)
/* Test engine API. */
#include "testeng_matop.h"

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extra, api, dimNum, align, loadFxn, procFxn ) { (fmt),extra,api,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }
#if 0 // HiFi3/3z API
static const tMatmulApi apimtx_mpy24x24      ={&mtx_mpy24x24_getScratchSize};
static const tMatmulApi apimtx_mpy24x24_fast ={&mtx_mpy24x24_fast_getScratchSize};
#endif
static const tMatmulApi apimtx_mpy16x16      ={&mtx_mpy16x16_getScratchSize};
static const tMatmulApi apimtx_mpyt16x16     ={&mtx_mpyt16x16_getScratchSize};
static const tMatmulApi apimtx_mpy32x32      ={&mtx_mpy32x32_getScratchSize};
static const tMatmulApi apimtx_mpyt32x32     ={&mtx_mpyt32x32_getScratchSize};
static const tMatmulApi apimtx_mpyf          ={&mtx_mpyf_getScratchSize};
static const tMatmulApi apimtx_mpy16x16_fast ={&mtx_mpy16x16_fast_getScratchSize};
static const tMatmulApi apimtx_mpy32x32_fast ={&mtx_mpy32x32_fast_getScratchSize};
static const tMatmulApi apimtx_mpyt16x16_fast={&mtx_mpyt16x16_fast_getScratchSize};
static const tMatmulApi apimtx_mpyt32x32_fast={&mtx_mpyt32x32_fast_getScratchSize};
static const tMatmulApi apimtx_mpyf_fast     ={&mtx_mpyf_fast_getScratchSize};
static const tMatmulApi apimtx_mpytf         ={&mtx_mpytf_getScratchSize};
static const tMatmulApi apimtx_mpytf_fast    ={&mtx_mpytf_fast_getScratchSize};
static const tMatmulApi apimtx_mpy8x8        ={&mtx_mpy8x8_getScratchSize};
static const tMatmulApi apimtx_mpyt8x8       ={&mtx_mpyt8x8_getScratchSize};
static const tMatmulApi apimtx_mpy8x8_fast   ={&mtx_mpy8x8_fast_getScratchSize};
static const tMatmulApi apimtx_mpyt8x8_fast  ={&mtx_mpyt8x8_fast_getScratchSize};
static const tMatmulApi apimtx_mpy8x16        ={&mtx_mpy8x16_getScratchSize};
static const tMatmulApi apimtx_mpyt8x16       ={&mtx_mpyt8x16_getScratchSize};
static const tMatmulApi apimtx_mpy8x16_fast   ={&mtx_mpy8x16_fast_getScratchSize};
static const tMatmulApi apimtx_mpyt8x16_fast  ={&mtx_mpyt8x16_fast_getScratchSize};

/* vec API test definitions. */
static const struct 
{
tTestEngTarget   funcList[MAX_FUNC_NUM];
tTestEngDesc     testDesc;
}
testDefTbl[] =
{
    { FUNC_LIST( (tTestEngTarget)&mtx_transpose16x16         ),TEST_DESC( FMT_REAL|FMT_FRACT16, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_NO,  &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transpose32x32         ),TEST_DESC( FMT_REAL|FMT_FRACT32, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_NO,  &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transpose8x8           ),TEST_DESC( FMT_REAL|FMT_INT8   , 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_NO,  &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transpose16x16_fast    ),TEST_DESC( FMT_REAL|FMT_FRACT16, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transpose32x32_fast    ),TEST_DESC( FMT_REAL|FMT_FRACT32, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transpose8x8_fast      ),TEST_DESC( FMT_REAL|FMT_INT8   , 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },

    /*
     * Stage 1
     */
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy16x16         ),TEST_DESC( FMT_REAL|FMT_FRACT16, MTX_PLAIN_LSH,&apimtx_mpy16x16      ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt16x16        ),TEST_DESC( FMT_REAL|FMT_FRACT16, MTX_PLAIN_LSH,&apimtx_mpyt16x16     ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy16x16      ),TEST_DESC( FMT_REAL|FMT_FRACT16, MTX_PLAIN_LSH,NULL                  ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy16x16_fast    ),TEST_DESC( FMT_REAL|FMT_FRACT16, MTX_PLAIN_LSH,&apimtx_mpy16x16_fast ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt16x16_fast   ),TEST_DESC( FMT_REAL|FMT_FRACT16, MTX_PLAIN_LSH,&apimtx_mpyt16x16_fast,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy16x16_fast ),TEST_DESC( FMT_REAL|FMT_FRACT16, MTX_PLAIN_LSH,NULL,                  TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
#if 0 // HiFi3/3z API
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy24x24 ),        TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,&apimtx_mpy24x24,     TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy24x24 ),     TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,             NULL,    TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy24x24_fast),    TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,&apimtx_mpy24x24_fast,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy24x24_fast ),TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,             NULL,    TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
#endif
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy32x32 ),         TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,&apimtx_mpy32x32      ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt32x32 ),        TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,&apimtx_mpyt32x32     ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy32x32 ),      TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,                  NULL,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy32x32_fast),     TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,&apimtx_mpy32x32_fast ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt32x32_fast),    TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,&apimtx_mpyt32x32_fast,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy32x32_fast ), TEST_DESC( FMT_REAL|FMT_FRACT32, MTX_PLAIN_LSH,                  NULL,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
 
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy8x8         ),TEST_DESC( FMT_REAL|FMT_INT8, MTX_PLAIN_LSH,&apimtx_mpy8x8      ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt8x8        ),TEST_DESC( FMT_REAL|FMT_INT8, MTX_PLAIN_LSH,&apimtx_mpyt8x8     ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy8x8      ),TEST_DESC( FMT_REAL|FMT_INT8, MTX_PLAIN_LSH,NULL                  ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy8x8_fast    ),TEST_DESC( FMT_REAL|FMT_INT8, MTX_PLAIN_LSH,&apimtx_mpy8x8_fast ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt8x8_fast   ),TEST_DESC( FMT_REAL|FMT_INT8, MTX_PLAIN_LSH,&apimtx_mpyt8x8_fast,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy8x8_fast ),TEST_DESC( FMT_REAL|FMT_INT8, MTX_PLAIN_LSH,NULL,                  TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
 
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy8x16         ),TEST_DESC( FMT_REAL|FMT_INT16, MTX_8X16_LSH,&apimtx_mpy8x16      ,TE_DIM_NUM_4, TE_ALIGN_NO,   &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt8x16        ),TEST_DESC( FMT_REAL|FMT_INT16, MTX_8X16_LSH,&apimtx_mpyt8x16     ,TE_DIM_NUM_4, TE_ALIGN_NO,   &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy8x16      ),TEST_DESC( FMT_REAL|FMT_INT16, MTX_8X16_LSH,NULL                  ,TE_DIM_NUM_4, TE_ALIGN_NO,  &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpy8x16_fast    ),TEST_DESC( FMT_REAL|FMT_INT16, MTX_8X16_LSH,&apimtx_mpy8x16_fast ,TE_DIM_NUM_4, TE_ALIGN_YES,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyt8x16_fast   ),TEST_DESC( FMT_REAL|FMT_INT16, MTX_8X16_LSH,&apimtx_mpyt8x16_fast,TE_DIM_NUM_4, TE_ALIGN_YES,  &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpy8x16_fast ),TEST_DESC( FMT_REAL|FMT_INT16, MTX_8X16_LSH,NULL,                  TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    /*
     * Stage 2
     */
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyf ),          TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpyf      ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpytf ),         TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpytf     ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpyf ),       TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL              ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyf_fast),      TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpyf_fast ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpytf_fast),     TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpytf_fast,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpyf_fast ),  TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL              ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transposef     ),TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_NO,  &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transposef_fast),TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },

    { FUNC_LIST( NULL ), TEST_DESC(  0,0, 0,0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for mat API functions. */
int main_matop( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int res = 1;

  #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), (seqFile),         \
                                                       isFull, isVerbose, breakOnError ) )
/*
   * Stage 1
   */
  if ( phaseNum == 0 || phaseNum == 1 )
  {
      DO_TEST( &mtx_mpy8x8        , "mtx_mpy8x8.seq"         );
      DO_TEST( &mtx_mpy8x8_fast   , "mtx_mpy8x8_fast.seq"    );
      DO_TEST( &mtx_mpyt8x8       , "mtx_mpyt8x8.seq"        );
      DO_TEST( &mtx_mpyt8x8_fast  , "mtx_mpyt8x8_fast.seq"   );
      DO_TEST( &mtx_vecmpy8x8     , "mtx_vecmpy8x8.seq"      );
      DO_TEST( &mtx_vecmpy8x8_fast, "mtx_vecmpy8x8_fast.seq" );

      DO_TEST( &mtx_mpy8x16        , "mtx_mpy8x16.seq"         );
      DO_TEST( &mtx_mpy8x16_fast   , "mtx_mpy8x16_fast.seq"    );
      DO_TEST( &mtx_mpyt8x16       , "mtx_mpyt8x16.seq"        );
      DO_TEST( &mtx_mpyt8x16_fast  , "mtx_mpyt8x16_fast.seq"   );
      DO_TEST( &mtx_vecmpy8x16     , "mtx_vecmpy8x16.seq"      );
      DO_TEST( &mtx_vecmpy8x16_fast, "mtx_vecmpy8x16_fast.seq" );

      DO_TEST( &mtx_mpy16x16        , "mtx_mpy16x16.seq"         );
      DO_TEST( &mtx_mpy16x16_fast   , "mtx_mpy16x16_fast.seq"    );
      DO_TEST( &mtx_mpyt16x16       , "mtx_mpyt16x16.seq"        );
      DO_TEST( &mtx_mpyt16x16_fast  , "mtx_mpyt16x16_fast.seq"   );
      DO_TEST( &mtx_vecmpy16x16     , "mtx_vecmpy16x16.seq"      );
      DO_TEST( &mtx_vecmpy16x16_fast, "mtx_vecmpy16x16_fast.seq" );
#if 0 // HiFi3/3z API
      DO_TEST( &mtx_mpy24x24        , "mtx_mpy24x24.seq"         );
      DO_TEST( &mtx_mpy24x24_fast   , "mtx_mpy24x24_fast.seq"    );
      DO_TEST( &mtx_vecmpy24x24     , "mtx_vecmpy24x24.seq"      );
      DO_TEST( &mtx_vecmpy24x24_fast, "mtx_vecmpy24x24_fast.seq" );
#endif
      DO_TEST( &mtx_mpy32x32        , "mtx_mpy32x32.seq"         );
      DO_TEST( &mtx_mpy32x32_fast   , "mtx_mpy32x32_fast.seq"    );
      DO_TEST( &mtx_mpyt32x32       , "mtx_mpyt32x32.seq"        );
      DO_TEST( &mtx_mpyt32x32_fast  , "mtx_mpyt32x32_fast.seq"   );
      DO_TEST( &mtx_vecmpy32x32     , "mtx_vecmpy32x32.seq"      );
      DO_TEST( &mtx_vecmpy32x32_fast, "mtx_vecmpy32x32_fast.seq" );

      DO_TEST( &mtx_transpose8x8          , "mtx_transpose8x8.seq"           );
      DO_TEST( &mtx_transpose8x8_fast     , "mtx_transpose8x8_fast.seq"      );
      DO_TEST( &mtx_transpose16x16        , "mtx_transpose16x16.seq"         );
      DO_TEST( &mtx_transpose16x16_fast   , "mtx_transpose16x16_fast.seq"    );
      DO_TEST( &mtx_transpose32x32        , "mtx_transpose32x32.seq"         );
      DO_TEST( &mtx_transpose32x32_fast   , "mtx_transpose32x32_fast.seq"    );
  }

  /*
   * Stage 2
   */

  if ( phaseNum == 0 || phaseNum == 2 )
  {
      DO_TEST( &mtx_mpyf             , "mtx_mpyf.seq"           );
      DO_TEST( &mtx_mpyf_fast        , "mtx_mpyf_fast.seq"      );
      DO_TEST( &mtx_mpytf            , "mtx_mpytf.seq"          );
      DO_TEST( &mtx_mpytf_fast       , "mtx_mpytf_fast.seq"     );
      DO_TEST( &mtx_vecmpyf          , "mtx_vecmpyf.seq"        );
      DO_TEST( &mtx_vecmpyf_fast     , "mtx_vecmpyf_fast.seq"   );
      DO_TEST( &mtx_transposef       , "mtx_transposef.seq"     );
      DO_TEST( &mtx_transposef_fast  , "mtx_transposef_fast.seq");
  }

  return (res);

} /* main_matop() */
