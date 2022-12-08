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
 * Test procedures for vector mathematics
 */

#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environemnt convfiguration. */
#include "config.h"
/* DSP Library API: arithmetic and logic functions on data vectors. */
#include LIBRARY_HEADER(math)

/* Test engine API. */
#include "testeng_math.h"


/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );


/* table for checking the bitexactness */
static const tTestEngVecSclTbl bitexactnessTbl[]=
{
#if 0 // HiFi3/3z API
    { te_math_processFxn_scl_recip32x32  , (te_fun_ptr_t)&scl_recip24x24      , (te_fun_ptr_t)&vec_recip24x24      },
    { te_math_processFxn_scl_divide32x32 , (te_fun_ptr_t)&scl_divide24x24     , (te_fun_ptr_t)&vec_divide24x24     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log2_24x24      , (te_fun_ptr_t)&vec_log2_24x24      },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_logn_24x24      , (te_fun_ptr_t)&vec_logn_24x24      },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log10_24x24     , (te_fun_ptr_t)&vec_log10_24x24     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog10_24x24 , (te_fun_ptr_t)&vec_antilog10_24x24 },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog2_24x24  , (te_fun_ptr_t)&vec_antilog2_24x24  },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilogn_24x24  , (te_fun_ptr_t)&vec_antilogn_24x24  },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sqrt24x24       , (te_fun_ptr_t)&vec_sqrt24x24       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sine24x24       , (te_fun_ptr_t)&vec_sine24x24       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_cosine24x24     , (te_fun_ptr_t)&vec_cosine24x24     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_atan24x24       , (te_fun_ptr_t)&vec_atan24x24       },
    { te_math_processFxn_scl_atan2       , (te_fun_ptr_t)&scl_atan2_24x24     , (te_fun_ptr_t)&vec_atan2_24x24     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tan24x24        , (te_fun_ptr_t)&vec_tan24x24        },
#endif
	// ( phaseNum == 0 || phaseNum == 1 )
    { te_math_processFxn_scl_recip16x16  , (te_fun_ptr_t)&scl_recip16x16      , (te_fun_ptr_t)&vec_recip16x16      },
    { te_math_processFxn_scl_recip32x32  , (te_fun_ptr_t)&scl_recip32x32      , (te_fun_ptr_t)&vec_recip32x32      },
    { te_math_processFxn_scl_divide16x16 , (te_fun_ptr_t)&scl_divide16x16     , (te_fun_ptr_t)&vec_divide16x16     },
    { te_math_processFxn_scl_divide32x32 , (te_fun_ptr_t)&scl_divide32x32     , (te_fun_ptr_t)&vec_divide32x32     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log2_32x32      , (te_fun_ptr_t)&vec_log2_32x32      },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_logn_32x32      , (te_fun_ptr_t)&vec_logn_32x32      },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log10_32x32     , (te_fun_ptr_t)&vec_log10_32x32     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog2_32x32  , (te_fun_ptr_t)&vec_antilog2_32x32  },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilogn_32x32  , (te_fun_ptr_t)&vec_antilogn_32x32  },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog10_32x32 , (te_fun_ptr_t)&vec_antilog10_32x32 },
    { te_math_processFxn_scl_vXvZ16      , (te_fun_ptr_t)&scl_sqrt32x16       , (te_fun_ptr_t)&vec_sqrt32x16       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sqrt32x32       , (te_fun_ptr_t)&vec_sqrt32x32       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sine32x32       , (te_fun_ptr_t)&vec_sine32x32       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_cosine32x32     , (te_fun_ptr_t)&vec_cosine32x32     },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_atan32x32       , (te_fun_ptr_t)&vec_atan32x32       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tan32x32        , (te_fun_ptr_t)&vec_tan32x32        },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sqrt16x16       , (te_fun_ptr_t)&vec_sqrt16x16       },
    { te_math_processFxn_scl_vXcvZ       , (te_fun_ptr_t)&scl_sqrt64x32       , (te_fun_ptr_t)&vec_sqrt64x32       },
    { te_math_processFxn_scl_recip16x16  , (te_fun_ptr_t)&scl_rsqrt16x16      , (te_fun_ptr_t)&vec_rsqrt16x16      },
    { te_math_processFxn_scl_recip32x32  , (te_fun_ptr_t)&scl_rsqrt32x32      , (te_fun_ptr_t)&vec_rsqrt32x32      },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tanh32x32       , (te_fun_ptr_t)&vec_tanh32x32       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sigmoid32x32    , (te_fun_ptr_t)&vec_sigmoid32x32    },
// ( phaseNum == 0 || phaseNum == 2 )
    { te_math_processFxn_scl_vX32sY32vZ  , (te_fun_ptr_t)&scl_int2float       , (te_fun_ptr_t)&vec_int2float       },
    { te_math_processFxn_scl_vXsY32vZ32  , (te_fun_ptr_t)&scl_float2int       , (te_fun_ptr_t)&vec_float2int       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sinef           , (te_fun_ptr_t)&vec_sinef           },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_cosinef         , (te_fun_ptr_t)&vec_cosinef         },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tanf            , (te_fun_ptr_t)&vec_tanf            },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log2f           , (te_fun_ptr_t)&vec_log2f           },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_lognf           , (te_fun_ptr_t)&vec_lognf           },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log10f          , (te_fun_ptr_t)&vec_log10f          },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog2f       , (te_fun_ptr_t)&vec_antilog2f       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilognf       , (te_fun_ptr_t)&vec_antilognf       },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog10f      , (te_fun_ptr_t)&vec_antilog10f      },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_atanf           , (te_fun_ptr_t)&vec_atanf           },
    { te_math_processFxn_scl_atan2       , (te_fun_ptr_t)&scl_atan2f          , (te_fun_ptr_t)&vec_atan2f          },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tanhf           , (te_fun_ptr_t)&vec_tanhf           },
    { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sigmoidf        , (te_fun_ptr_t)&vec_sigmoidf    },
    { te_math_processFxn_scl_vXvYvZ      , (te_fun_ptr_t)&scl_powf            , (te_fun_ptr_t)&vec_powf        },

	{0, NULL, NULL},
};


#define DO_TEST( fxn, seqFile, extraFlags )                                                         \
if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), (seqFile),             \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) )
/* Extended variant runs each SEQ-file twice, with extended Error Handling check in the second
* invocation. */
#define DO_TEST_EXT( fxn, seqFile, extraFlags )                                                   \
    { if ( res || !breakOnError ) res &= ( 0 != testExec(                                         \
                                                    (tTestEngTarget)(fxn),                        \
                                                    (seqFile),                                    \
                                                    (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                    isFull, isVerbose, breakOnError, 0 ) );       \
    if ( res || !breakOnError ) res &= ( 0 != testExec(                                           \
                                                    (tTestEngTarget)(fxn),                        \
                                                    (seqFile),                                    \
                                                    (extraFlags) | TE_ERRH_EXTENDED_TEST_ENABLE,  \
                                                    isFull, isVerbose, breakOnError, 0 ) ); }

#define DO_TEST_BITEXACTNESS(fxn, seqFile, extraFlags )                                               \
{                                                                                                     \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), (seqFile),             \
                                            (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,             \
                                         isFull, isVerbose, breakOnError,0 ) );                       \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), (seqFile),             \
                                            (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,             \
                                         isFull, isVerbose, breakOnError,1 ) ); }

/* Perform all tests for math API functions. */
int main_mathv( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int res = 1;
/*
   * Stage 1
   */
  if ( phaseNum == 0 || phaseNum == 1 )
  {
        DO_TEST             ( &vec_divide16x16        , "vec_divide16x16.seq"      , 0 );
#if 0 //HiFi3/3z API
        DO_TEST             ( &vec_divide24x24        , "vec_divide24x24.seq"      , 0 );
#endif
        DO_TEST             ( &vec_divide32x32        , "vec_divide32x32.seq"      , 0 );
        DO_TEST             ( &vec_divide64x32i       , "vec_divide64x32i.seq"     , 0 );
        DO_TEST             ( &vec_divide64x64        , "vec_divide64x64.seq"      , 0 );
        DO_TEST             ( &vec_recip16x16         , "vec_recip16x16.seq"       , 0 );
#if 0 //HiFi3/3z API
        DO_TEST             ( &vec_recip24x24         , "vec_recip24x24.seq"       , 0 );
#endif
        DO_TEST             ( &vec_recip32x32         , "vec_recip32x32.seq"       , 0 );
        DO_TEST             ( &vec_recip64x64         , "vec_recip64x64.seq"       , 0 );
        DO_TEST_BITEXACTNESS( &vec_log2_32x32         , "vec_log2_32x32.seq"       , 0 );
        DO_TEST_BITEXACTNESS( &vec_logn_32x32         , "vec_logn_32x32.seq"       , 0 );
        DO_TEST_BITEXACTNESS( &vec_log10_32x32        , "vec_log10_32x32.seq"      , 0 );
#if 0 //HiFi3/3z API
        DO_TEST_BITEXACTNESS( &vec_log2_24x24         , "vec_log2_24x24.seq"       , 0 );
        DO_TEST_BITEXACTNESS( &vec_logn_24x24         , "vec_logn_24x24.seq"       , 0 );
        DO_TEST_BITEXACTNESS( &vec_log10_24x24        , "vec_log10_24x24.seq"      , 0 );
#endif
        DO_TEST_BITEXACTNESS( &vec_antilog2_32x32     , "vec_antilog2_32x32.seq"   , 0 );
        DO_TEST_BITEXACTNESS( &vec_antilogn_32x32     , "vec_antilogn_32x32.seq"   , 0 );
        DO_TEST_BITEXACTNESS( &vec_antilog10_32x32    , "vec_antilog10_32x32.seq"  , 0 );
#if 0 //HiFi3/3z API
        DO_TEST_BITEXACTNESS( &vec_antilog10_24x24    , "vec_antilog10_24x24.seq"  , 0 );
        DO_TEST_BITEXACTNESS( &vec_antilog2_24x24     , "vec_antilog2_24x24.seq"   , 0 );
        DO_TEST_BITEXACTNESS( &vec_antilogn_24x24     , "vec_antilogn_24x24.seq"   , 0 );
#endif
        DO_TEST             ( &vec_pow_32x32          , "vec_pow_32x32.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_sqrt32x32          , "vec_sqrt32x32.seq"        , 0 );
#if 0 //HiFi3/3z API
        DO_TEST_BITEXACTNESS( &vec_sqrt24x24          , "vec_sqrt24x24.seq"        , 0 );
#endif
        DO_TEST_BITEXACTNESS( &vec_sqrt16x16          , "vec_sqrt16x16.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_sqrt32x16          , "vec_sqrt32x16.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_sqrt64x32          , "vec_sqrt64x32.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_sine32x32          , "vec_sine32x32.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_cosine32x32        , "vec_cosine32x32.seq"      , 0 );
#if 0 //HiFi3/3z API
        DO_TEST_BITEXACTNESS( &vec_sine24x24          , "vec_sine24x24.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_cosine24x24        , "vec_cosine24x24.seq"      , 0 );
#endif
        DO_TEST_BITEXACTNESS( &vec_atan32x32          , "vec_atan32x32.seq"        , 0 );
#if 0 //HiFi3/3z API
        DO_TEST_BITEXACTNESS( &vec_atan24x24          , "vec_atan24x24.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_atan2_24x24        , "vec_atan2_24x24.seq"      , 0 );
#endif
        DO_TEST_BITEXACTNESS( &vec_tan32x32           , "vec_tan32x32.seq"         , 0 );
#if 0 //HiFi3/3z API
        DO_TEST_BITEXACTNESS( &vec_tan24x24           , "vec_tan24x24.seq"         , 0 );
#endif
        DO_TEST             ( &vec_rsqrt16x16         , "vec_rsqrt16x16.seq"       , 0 );
        DO_TEST             ( &vec_rsqrt32x32         , "vec_rsqrt32x32.seq"       , 0 );
        DO_TEST_BITEXACTNESS( &vec_tanh32x32          , "vec_tanh32x32.seq"        , 0 );
        DO_TEST_BITEXACTNESS( &vec_sigmoid32x32       , "vec_sigmoid32x32.seq"     , 0 );
        DO_TEST             ( &vec_softmax32x32       , "vec_softmax32x32.seq"     , 0 );
        DO_TEST             ( &vec_relu32x32          , "vec_relu32x32.seq"        , 0 );
  }

  /*
   * Stage 2
   */
  if ( phaseNum == 0 || phaseNum == 2 )
  {
       DO_TEST_BITEXACTNESS( &vec_int2float     , "vec_int2float.seq"       , 0 );
       DO_TEST_BITEXACTNESS( &vec_float2int     , "vec_float2int.seq"       , 0 );
       DO_TEST_BITEXACTNESS( &vec_sinef         , "vec_sinef.seq"           , 0 );
       DO_TEST_BITEXACTNESS( &vec_cosinef       , "vec_cosinef.seq"         , 0 );
       DO_TEST_BITEXACTNESS( &vec_tanf          , "vec_tanf.seq"            , 0 );
       DO_TEST_BITEXACTNESS( &vec_log2f         , "vec_log2f.seq"           , 0 );
       DO_TEST_BITEXACTNESS( &vec_lognf         , "vec_lognf.seq"           , 0 );
       DO_TEST_BITEXACTNESS( &vec_log10f        , "vec_log10f.seq"          , 0 );
       DO_TEST_BITEXACTNESS( &vec_antilog2f     , "vec_antilog2f.seq"       , 0 );
       DO_TEST_BITEXACTNESS( &vec_antilognf     , "vec_antilognf.seq"       , 0 );
       DO_TEST_BITEXACTNESS( &vec_antilog10f    , "vec_antilog10f.seq"      , 0 );
       DO_TEST_BITEXACTNESS( &vec_powf          , "vec_powf.seq"            , 0 );
       DO_TEST_BITEXACTNESS( &vec_atanf         , "vec_atanf.seq"           , 0 );
       DO_TEST_BITEXACTNESS( &vec_atan2f        , "vec_atan2f.seq"          , 0 );
       DO_TEST_BITEXACTNESS( &vec_tanhf         , "vec_tanhf.seq"           , 0 );
       DO_TEST_BITEXACTNESS( &vec_sigmoidf      , "vec_sigmoidf.seq"        , 0 );
       DO_TEST             ( &vec_softmaxf      , "vec_softmaxf.seq"        , 0 );
       DO_TEST             ( &vec_reluf         , "vec_reluf.seq"           , 0 );

  }

  return (res);

} /* main_mathv() */

/* Perform all tests for math API functions. */
int main_mathvf( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int res = 1;

  /*
   * Stage 1
   */
  if ( phaseNum == 0 || phaseNum == 1 )
  {
        DO_TEST( &vec_divide16x16_fast   , "vec_divide16x16_fast.seq" , 0 );
        DO_TEST( &vec_divide32x32_fast   , "vec_divide32x32_fast.seq" , 0 );
#if 0 //HiFi3/3z API
        DO_TEST( &vec_divide24x24_fast   , "vec_divide24x24_fast.seq" , 0 );
#endif
        DO_TEST( &vec_sqrt32x32_fast     , "vec_sqrt32x32_fast.seq"   , 0 );
#if 0 //HiFi3/3z API
        DO_TEST( &vec_sqrt24x24_fast     , "vec_sqrt24x24_fast.seq"   , 0 );
#endif
        DO_TEST( &vec_sine32x32_fast     , "vec_sine32x32_fast.seq"   , 0 );
        DO_TEST( &vec_cosine32x32_fast   , "vec_cosine32x32_fast.seq" , 0 );
#if 0 //HiFi3/3z API
        DO_TEST( &vec_sine24x24_fast     , "vec_sine24x24_fast.seq"   , 0 );
        DO_TEST( &vec_cosine24x24_fast   , "vec_cosine24x24_fast.seq" , 0 );
#endif
  }

  /*
   * Stage 2
   */
  if ( phaseNum == 0 || phaseNum == 2 )
  {
  }
  return (res);
} /* main_mathvf() */

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget   targetFxn, const char * seqName, 
              int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness )
{
  #define MAX_FUNC_NUM   16
  /* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
  #define FUNC_LIST(...) { __VA_ARGS__, NULL }
  /* Initializer for a test description structure. */
  #define TEST_DESC( fmt,extraParam,dimNum, align, loadFxn, procFxn ) { (fmt),extraParam,bitexactnessTbl,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

  /* vec API test definitions. */
  static const struct 
  {
    tTestEngTarget   funcList[MAX_FUNC_NUM];
    tTestEngDesc     testDesc;
  }
  testDefTbl[] =
  {
    /*
     * Stage 1
     */
#if 0 //HiFi3/3z API
    {
      FUNC_LIST( (tTestEngTarget)&vec_atan24x24),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sine24x24_fast,(tTestEngTarget)&vec_cosine24x24_fast),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sine24x24,     (tTestEngTarget)&vec_cosine24x24,   (tTestEngTarget)&vec_tan24x24,
                 (tTestEngTarget)&vec_log2_24x24,    (tTestEngTarget)&vec_logn_24x24,    (tTestEngTarget)&vec_log10_24x24,
                 (tTestEngTarget)&vec_antilog2_24x24,(tTestEngTarget)&vec_antilogn_24x24,(tTestEngTarget)&vec_antilog10_24x24),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( ((tTestEngTarget)&vec_sqrt24x24_fast),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( ((tTestEngTarget)&vec_sqrt24x24),
      TEST_DESC( FMT_REAL|FMT_FRACT32,OVLP_XZ,  TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_recip24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZd, &te_math_processFxn_recip32x32) },
    {
      FUNC_LIST( ((tTestEngTarget)&vec_divide24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, OVLP_XZ | OVLP_YZ ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZd, &te_math_processFxn_divide32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide24x24_fast),
      TEST_DESC(FMT_REAL | FMT_FRACT32, 0,TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZd, &te_math_processFxn_divide32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_atan2_24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, OVLP_XZ|OVLP_YZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvYvX) },
#endif
    {
      FUNC_LIST( (tTestEngTarget)&vec_sine32x32,     (tTestEngTarget)&vec_cosine32x32,   (tTestEngTarget)&vec_tan32x32,
                 (tTestEngTarget)&vec_log2_32x32,    (tTestEngTarget)&vec_logn_32x32,    (tTestEngTarget)&vec_log10_32x32,
                 (tTestEngTarget)&vec_antilog2_32x32,(tTestEngTarget)&vec_antilogn_32x32,(tTestEngTarget)&vec_antilog10_32x32),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_atan32x32),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sine32x32_fast,(tTestEngTarget)&vec_cosine32x32_fast),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sqrt32x32),
      TEST_DESC( FMT_REAL|FMT_FRACT32,OVLP_XZ,  TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sqrt32x32_fast),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_recip16x16),
      TEST_DESC( FMT_REAL|FMT_FRACT16, OVLP_XZ | OVLP_XW,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZf, &te_math_processFxn_recip16x16 ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_recip32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZd, &te_math_processFxn_recip32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide16x16),
      TEST_DESC(FMT_REAL | FMT_FRACT16, OVLP_XZ | OVLP_XW | OVLP_YZ | OVLP_YW,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZf, &te_math_processFxn_divide16x16) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide16x16_fast),
      TEST_DESC(FMT_REAL | FMT_FRACT16, OVLP_XZ | OVLP_XW | OVLP_YZ | OVLP_YW,TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZf, &te_math_processFxn_divide16x16) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, OVLP_XZ | OVLP_YZ ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZd, &te_math_processFxn_divide32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide32x32_fast),
      TEST_DESC(FMT_REAL | FMT_FRACT32, 0,TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZd, &te_math_processFxn_divide32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide64x32i),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_YZ, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vX64vYvZ, &te_processFxn_vZvXvY ) },

    {
      FUNC_LIST( (tTestEngTarget)&vec_rsqrt32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, OVLP_XZ ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZd, &te_math_processFxn_recip32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_rsqrt16x16),
      TEST_DESC(FMT_REAL | FMT_FRACT16, OVLP_XZ | OVLP_XW,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZf, &te_math_processFxn_recip16x16) },

    {
      FUNC_LIST( (tTestEngTarget)&vec_tanh32x32, (tTestEngTarget)&vec_sigmoid32x32, (tTestEngTarget)&vec_softmax32x32),
      TEST_DESC( FMT_FRACT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sqrt16x16),
      TEST_DESC( FMT_FRACT16,OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sqrt64x32),
      TEST_DESC( FMT_FRACT32,0, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vX64vZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_relu32x32 ),
      TEST_DESC( FMT_FRACT32, 0*OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZ_vXsY ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sqrt32x16),
      TEST_DESC( FMT_REAL|FMT_FRACT32,0,  TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vX32vZ16, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_recip64x64),
      TEST_DESC(FMT_REAL | FMT_INT64, 0,TE_DIM_NUM_1, TE_ALIGN_NO, &te_math_loadFxn_recip64x64, &te_math_processFxn_vec_recip64x64) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide64x64),
      TEST_DESC(FMT_REAL | FMT_INT64, 0,TE_DIM_NUM_1, TE_ALIGN_NO, &te_math_loadFxn_divide64x64, &te_math_processFxn_vec_divide64x64) },

    {
      FUNC_LIST( (tTestEngTarget)&vec_pow_32x32),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ|OVLP_YZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZdp, &te_math_processFxn_pow32x32) },
    /*
     * Stage 2
     */
    {
        FUNC_LIST((tTestEngTarget)&vec_powf),
        TEST_DESC(FMT_REAL | FMT_FLOAT32,OVLP_XZ|OVLP_YZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY) },
    { 
      FUNC_LIST( (tTestEngTarget)&vec_int2float ),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vX32sY32vZ, &te_processFxn_vZvXsY32 ) },
    { 
      FUNC_LIST( (tTestEngTarget)&vec_float2int ),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32vZ32, &te_processFxn_vZvXsY32 ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sinef,(tTestEngTarget)&vec_cosinef,(tTestEngTarget)&vec_tanf,
                 (tTestEngTarget)&vec_log2f,(tTestEngTarget)&vec_lognf,(tTestEngTarget)&vec_log10f,
                 (tTestEngTarget)&vec_atanf,(tTestEngTarget)&vec_antilog2f,(tTestEngTarget)&vec_antilognf,
                 (tTestEngTarget)&vec_antilog10f),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_atan2f ),
      TEST_DESC( FMT_FLOAT32, OVLP_XZ|OVLP_YZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvYvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_tanhf, (tTestEngTarget)&vec_sigmoidf, (tTestEngTarget)&vec_softmaxf),
      TEST_DESC( FMT_FLOAT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_reluf ),
      TEST_DESC( FMT_FLOAT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZ_vXsY ) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, 0, NULL, NULL ) } /* End of table */
  };

  {
    int tblIx, funcIx;

    for ( tblIx=0; tblIx<(int)(sizeof(testDefTbl)/sizeof(testDefTbl[0])); tblIx++ )
    {
      for ( funcIx=0; funcIx<MAX_FUNC_NUM; funcIx++ )
      {
        if ( targetFxn == testDefTbl[tblIx].funcList[funcIx] )
        {
          tTestEngDesc testDesc = testDefTbl[tblIx].testDesc;
          testDesc.extraParam |= (uint32_t)errhExtendedTest;

          return ( TestEngRun( targetFxn, &testDesc, 
                               seqName, isFull, 
                               isVerbose, breakOnError, testBitexactness ) );
        }
      }
    }

    ASSERT( !"Test not defined" );
    return (0);
  }
  return te_Exec(testDefTbl, sizeof(testDefTbl) / sizeof(testDefTbl[0]), MAX_FUNC_NUM, targetFxn, seqName, isFull, isVerbose, breakOnError);

} /* testExec() */
