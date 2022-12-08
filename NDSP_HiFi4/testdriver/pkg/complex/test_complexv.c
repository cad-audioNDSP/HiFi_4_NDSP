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
/* Library API */
#include LIBRARY_HEADER(complex)
/* Test engine API. */
#include "testeng.h"

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );

/* Apply a function to the test case data set:
 * scalar functions with single argument, e.g. cos(). Input X is complex, output is real*/
extern void processFxn_scl_vXcvZ( tTestEngContext * context );

/* table for checking the bitexactness */
static const tTestEngVecSclTbl bitexactnessTbl[]=
{
// ( phaseNum == 0 || phaseNum == 1 )

// ( phaseNum == 0 || phaseNum == 2 )
    { processFxn_scl_vXcvZ , (te_fun_ptr_t)&scl_complex2mag    , (te_fun_ptr_t)&vec_complex2mag    },
    { processFxn_scl_vXcvZ , (te_fun_ptr_t)&scl_complex2invmag , (te_fun_ptr_t)&vec_complex2invmag },

    {0, NULL, NULL},
};

int main_complexv( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int res = 1;
  #define DO_TEST( fxn, seqFile, extraFlags )                                                         \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), (seqFile),             \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) )
  /* Extended variant runs each SEQ-file twice, with extended Error Handling check in the second
   * invocation. */
  #define DO_TEST_EXT( fxn, seqFile, extraFlags )                                                     \
      { if ( res || !breakOnError ) res &= ( 0 != testExec(                                           \
                                                        (tTestEngTarget)(fxn),                        \
                                                        (seqFile),                                    \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) );       \
        if ( res || !breakOnError ) res &= ( 0 != testExec(                                           \
                                                        (tTestEngTarget)(fxn),                        \
                                                        (seqFile),                                    \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_ENABLE,  \
                                                        isFull, isVerbose, breakOnError, 0 ) ); }

#define DO_TEST_BITEXACTNESS(fxn, seqFile, extraFlags)                                                \
{                                                                                                     \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), (seqFile),             \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                         isFull, isVerbose, breakOnError,0 ) );                       \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), (seqFile),             \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                         isFull, isVerbose, breakOnError,1 ) ); }
  /*
   * Stage 1
   */
  if ( phaseNum == 0 || phaseNum == 1 )
  {
  }

  /*
   * Stage 2
   */
  if ( phaseNum == 0 || phaseNum == 2 )
  {
        DO_TEST_BITEXACTNESS ( &vec_complex2mag   , "vec_complex2mag.seq"     , 0 );
        DO_TEST_BITEXACTNESS ( &vec_complex2invmag, "vec_complex2invmag.seq"  , 0 );
  }

  return (res);

} /* main_complexv() */

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
    /*
     * Stage 2
     */

    {
      FUNC_LIST( (tTestEngTarget)&vec_complex2mag,(tTestEngTarget)&vec_complex2invmag ),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXcvZ, &te_processFxn_vZvX ) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0,0, 0, NULL, NULL ) } /* End of table */
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
          testDesc.extraParam = (uint32_t)errhExtendedTest;

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
