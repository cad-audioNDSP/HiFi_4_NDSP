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
 * Test engine extension for FFT tests
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
/* Test engine extension for DCT */
#include "testeng_dct.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test environment utils. */
#include "utils.h"

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Suppress Visual C warnings on +/-INFINITY macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

#define sz_fp64c    sizeof(complex_double)

typedef int (*dct_fxns_t)( void* y,void* x,dct_handle_t h,int scalingOpt);
typedef int (*dct_fxns_2d_t)( void* y,void* x,dct_handle_t h,int L,int scalingOpt);
typedef int (*dct_fxnsf_t)( void* y,void* x,dct_handle_t h);

typedef struct
{
    int           N;
    dct_fxns_t    fun;
    dct_handle_t  h;
}
tHandleTbl;


/* find dct handle */
void findItem(tHandleTbl* tbl, int N, dct_fxns_t fun)
{
    const tHandleTbl handle_tbl[]=
    {
        // DCT-II
        { 32, (dct_fxns_t)dct_16x16, (dct_handle_t)dct2_16_32},
        { 64, (dct_fxns_t)dct_16x16, (dct_handle_t)dct2_16_64},
        { 32, (dct_fxns_t)dct_32x16, (dct_handle_t)dct2_16_32},
        { 64, (dct_fxns_t)dct_32x16, (dct_handle_t)dct2_16_64},
#if 0 // HiFi3/3z API
        { 32, (dct_fxns_t)dct_24x24, (dct_handle_t)dct2_32_32},
        { 64, (dct_fxns_t)dct_24x24, (dct_handle_t)dct2_32_64},
#else
        { 32, (dct_fxns_t)NULL, (dct_handle_t)dct2_32_32},
        { 64, (dct_fxns_t)NULL, (dct_handle_t)dct2_32_64},
#endif
        { 32, (dct_fxns_t)dct_32x32, (dct_handle_t)dct2_32_32},
        { 64, (dct_fxns_t)dct_32x32, (dct_handle_t)dct2_32_64},
        { 32, (dct_fxns_t)dctf     , (dct_handle_t)dct2_f_32 },
        { 64, (dct_fxns_t)dctf     , (dct_handle_t)dct2_f_64 },
    #if 1
        // DCT-IV handles
        { 32 , (dct_fxns_t)dct4_32x16, (dct_handle_t)dct4_16_32 },
        { 64 , (dct_fxns_t)dct4_32x16, (dct_handle_t)dct4_16_64 },
        { 128, (dct_fxns_t)dct4_32x16, (dct_handle_t)dct4_16_128},
        { 256, (dct_fxns_t)dct4_32x16, (dct_handle_t)dct4_16_256},
        { 512, (dct_fxns_t)dct4_32x16, (dct_handle_t)dct4_16_512},
        { 32 , (dct_fxns_t)dct4_32x32, (dct_handle_t)dct4_32_32 },
        { 64 , (dct_fxns_t)dct4_32x32, (dct_handle_t)dct4_32_64 },
        { 128, (dct_fxns_t)dct4_32x32, (dct_handle_t)dct4_32_128},
        { 256, (dct_fxns_t)dct4_32x32, (dct_handle_t)dct4_32_256},
        { 512, (dct_fxns_t)dct4_32x32, (dct_handle_t)dct4_32_512},
#if 0 // HiFi3/3z API
        { 32 , (dct_fxns_t)dct4_24x24, (dct_handle_t)dct4_32_32 },
        { 64 , (dct_fxns_t)dct4_24x24, (dct_handle_t)dct4_32_64 },
        { 128, (dct_fxns_t)dct4_24x24, (dct_handle_t)dct4_32_128},
        { 256, (dct_fxns_t)dct4_24x24, (dct_handle_t)dct4_32_256},
        { 512, (dct_fxns_t)dct4_24x24, (dct_handle_t)dct4_32_512},
#else
        { 32 , (dct_fxns_t)NULL, (dct_handle_t)dct4_32_32 },
        { 64 , (dct_fxns_t)NULL, (dct_handle_t)dct4_32_64 },
        { 128, (dct_fxns_t)NULL, (dct_handle_t)dct4_32_128},
        { 256, (dct_fxns_t)NULL, (dct_handle_t)dct4_32_256},
        { 512, (dct_fxns_t)NULL, (dct_handle_t)dct4_32_512},
#endif
    #endif
    #if 1
        // MDCT/IMDCT handles
        { 32 , (dct_fxns_t)mdct_32x16, (dct_handle_t)mdct_16_32 },
        { 64 , (dct_fxns_t)mdct_32x16, (dct_handle_t)mdct_16_64 },
        { 128, (dct_fxns_t)mdct_32x16, (dct_handle_t)mdct_16_128},
        { 256, (dct_fxns_t)mdct_32x16, (dct_handle_t)mdct_16_256},
        { 512, (dct_fxns_t)mdct_32x16, (dct_handle_t)mdct_16_512},
        { 32 , (dct_fxns_t)mdct_32x32, (dct_handle_t)mdct_32_32 },
        { 64 , (dct_fxns_t)mdct_32x32, (dct_handle_t)mdct_32_64 },
        { 128, (dct_fxns_t)mdct_32x32, (dct_handle_t)mdct_32_128},
        { 256, (dct_fxns_t)mdct_32x32, (dct_handle_t)mdct_32_256},
        { 512, (dct_fxns_t)mdct_32x32, (dct_handle_t)mdct_32_512},
#if 0 // HiFi3/3z API
        { 32 , (dct_fxns_t)mdct_24x24, (dct_handle_t)mdct_32_32 },
        { 64 , (dct_fxns_t)mdct_24x24, (dct_handle_t)mdct_32_64 },
        { 128, (dct_fxns_t)mdct_24x24, (dct_handle_t)mdct_32_128},
        { 256, (dct_fxns_t)mdct_24x24, (dct_handle_t)mdct_32_256},
        { 512, (dct_fxns_t)mdct_24x24, (dct_handle_t)mdct_32_512},
#else
        { 32 , (dct_fxns_t)NULL, (dct_handle_t)mdct_32_32 },
        { 64 , (dct_fxns_t)NULL, (dct_handle_t)mdct_32_64 },
        { 128, (dct_fxns_t)NULL, (dct_handle_t)mdct_32_128},
        { 256, (dct_fxns_t)NULL, (dct_handle_t)mdct_32_256},
        { 512, (dct_fxns_t)NULL, (dct_handle_t)mdct_32_512},
#endif
        { 32 , (dct_fxns_t)imdct_32x16, (dct_handle_t)mdct_16_32 },
        { 64 , (dct_fxns_t)imdct_32x16, (dct_handle_t)mdct_16_64 },
        { 128, (dct_fxns_t)imdct_32x16, (dct_handle_t)mdct_16_128},
        { 256, (dct_fxns_t)imdct_32x16, (dct_handle_t)mdct_16_256},
        { 512, (dct_fxns_t)imdct_32x16, (dct_handle_t)mdct_16_512},
        { 32 , (dct_fxns_t)imdct_32x32, (dct_handle_t)mdct_32_32 },
        { 64 , (dct_fxns_t)imdct_32x32, (dct_handle_t)mdct_32_64 },
        { 128, (dct_fxns_t)imdct_32x32, (dct_handle_t)mdct_32_128},
        { 256, (dct_fxns_t)imdct_32x32, (dct_handle_t)mdct_32_256},
        { 512, (dct_fxns_t)imdct_32x32, (dct_handle_t)mdct_32_512},
#if 0 // HiFi3/3z API
        { 32 , (dct_fxns_t)imdct_24x24, (dct_handle_t)mdct_32_32 },
        { 64 , (dct_fxns_t)imdct_24x24, (dct_handle_t)mdct_32_64 },
        { 128, (dct_fxns_t)imdct_24x24, (dct_handle_t)mdct_32_128},
        { 256, (dct_fxns_t)imdct_24x24, (dct_handle_t)mdct_32_256},
        { 512, (dct_fxns_t)imdct_24x24, (dct_handle_t)mdct_32_512},
#else
        { 32 , (dct_fxns_t)NULL, (dct_handle_t)mdct_32_32 },
        { 64 , (dct_fxns_t)NULL, (dct_handle_t)mdct_32_64 },
        { 128, (dct_fxns_t)NULL, (dct_handle_t)mdct_32_128},
        { 256, (dct_fxns_t)NULL, (dct_handle_t)mdct_32_256},
        { 512, (dct_fxns_t)NULL, (dct_handle_t)mdct_32_512},
#endif
    #endif
    #if 1
        // DCT-2D handles
        { 8 , (dct_fxns_t)dct2d_8x16, (dct_handle_t)dct2d_16_8 },
        { 8 , (dct_fxns_t)idct2d_16x8, (dct_handle_t)idct2d_16_8 },
    #endif
        {0,0,0}
    };
    int n;
    const tHandleTbl* pTbl=handle_tbl;
    for (n=0; n<(int)(sizeof(handle_tbl)/sizeof(handle_tbl[0])); n++,pTbl++)
    {
        if(pTbl->N==N && pTbl->fun==fun) 
        {
            tbl[0]=pTbl[0];
            return;
        }
    }
    tbl->fun=NULL;
    tbl->h=NULL;
    tbl->N=0;
    ASSERT(0);  //  not found!
}


/* FFT test context. */
typedef struct tagTestEngContext_dct_int
{
  tTestEngContext_dct ext; /* Externally visible part of FFT context. */

} tTestEngContext_dct_int;

/* Create a target algorithm instance and set tTestEngContext::target fields.
 * Return zero if failed. */
int te_create_dct( tTestEngContext * context )
{
  tTestEngContext_dct_int * context_dct;

  int res;

  /*
   * Allocate and initialize the context structure.
   */

  context_dct = (tTestEngContext_dct_int * )malloc( sizeof(*context_dct) );

  if ( !( res = ( 0 != context_dct ) ) )
  {
    printf( "te_create_dct(): malloc() failed\n" );
  }
  else
  {
    memset( context_dct, 0, sizeof(*context_dct) );
  }

  /*
   * Load twiddle factor tables.
   */
  if ( res )
  {

    context->target.handle = &context_dct->ext;
  }

  {
    typedef int (*tFxn)(void * y, const void * x, int N );
    tFxn fxn_fwd,fxn_inv ;
    fxn_fwd = (tFxn)((const tTestEngDesc_fft *)context->desc)->frwTransFxn;
    fxn_inv = (tFxn)((const tTestEngDesc_fft *)context->desc)->invTransFxn;
    if(!IS_PRESENT(fxn_fwd) &&
       !IS_PRESENT(fxn_inv))
    {
        // FUT is not defined
        return -1;
    }
  }

  return (res);

} /* te_create_dct() */

/* Destroy the target algorithm instance and free memory block(s) allocated
 * for the target object. Return zero whenever any violation of memory area
 * bounds is detected. */
int te_destroy_dct( tTestEngContext * context )
{
  tTestEngContext_dct_int * context_dct;

  ASSERT( context );

  if ( 0 != ( context_dct = (tTestEngContext_dct_int *)context->target.handle ) )
  {
    free( context_dct );
  }

  return (1);

} /* te_destroy_dct() */

/* Allocate in/out vectors for the next test case, and load the data set
 * from the SEQ-file. Return zero if failed. */
int te_load_dct( tTestEngContext * context )
{
  tTestEngContext_dct_int * context_dct = (tTestEngContext_dct_int *)context->target.handle;
  tVec BEXP, Z, Zlo, Zhi;
  int res = 0;
  int isFract = ( ( FMT_FRACT16 == ( context->desc->fmt & FMT_DTYPE_MASK ) ) ||
                  ( FMT_FRACT32 == ( context->desc->fmt & FMT_DTYPE_MASK ) ) );

  NASSERT( context_dct );

  memset( &context_dct->ext, 0, sizeof(context_dct->ext) );

  memset( &BEXP, 0, sizeof(BEXP) );
  memset( &Z   , 0, sizeof(Z   ) );
  memset( &Zlo , 0, sizeof(Zlo ) );
  memset( &Zhi , 0, sizeof(Zhi ) );

  /* If DCT supports the scaling option, read the scaling method from the SEQ-file. */
  if ( ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) &&
       ( 1 != seqFileScanf( context->seqFile, "%d", &context_dct->ext.scale_method ) ) )
  {
    printf( "te_load_dct(): bad SEQ-file format (a)\n" );
  }
  /* For a fixed point blockwise DCT, allocate a vector for temporal storage of block exponent. */
  else if ( 0 != ( context->desc->extraParam & TE_FFT_BLOCKWISE ) && isFract &&
            !vecAlloc( &BEXP, context->args.dim[2], TE_ALIGN_NO, FMT_INT16, 0 ) )
  {
    printf( "te_load_dct(): failed to allocate BEXP, L=%d\n", context->args.dim[2] );
  }
  /* Read input data filename. */
  else if ( 1 != seqFileScanf( context->seqFile, "%63s", &context_dct->ext.fInName ) )
  {
    printf( "te_load_dct(): bad SEQ-file format (b)\n" );
  }
  /* Read reference data filename. */
  else if ( 1 != seqFileScanf( context->seqFile, "%63s", &context_dct->ext.fRefName ) )
  {
    printf( "te_load_dct(): bad SEQ-file format (c)\n" );
  }
  /* Allocate vectors for SINAD verification. */
  else if ( 3 != vecsAlloc( TE_ALIGN_NO, FMT_FLOAT32, &Z, 1, &Zlo, 1, &Zhi, 1, 0 ) )
  {
    printf( "te_load_dct(): failed to allocate vectors Z/Zlo/Zhi\n" );
  }
  /* Read the minimum SINAD value from the SEQ-file. */
  else if ( 1 != seqFileScanf( context->seqFile, "%f", vecGetElem_fl32( &Zlo, 0 ) ) )
  {
    printf( "te_load_dct(): bad SEQ-file format (d)\n" );
  }
  else
  {
    /* Set SINAD upper limit to infinity. */
    *vecGetElem_fl32( &Zhi, 0 ) = INFINITY;

    memset( &context->dataSet, 0, sizeof(context->dataSet) );

    context->dataSet.X   = BEXP;
    context->dataSet.Z   = Z;
    context->dataSet.Zlo = Zlo;
    context->dataSet.Zhi = Zhi;

    res = 1;
  }

  if ( !res )
  {
    if ( BEXP.szBulk ) vecFree( &BEXP );
    if ( Z   .szBulk ) vecFree( &Z    );
    if ( Zlo .szBulk ) vecFree( &Zlo  );
    if ( Zhi .szBulk ) vecFree( &Zhi  );
  }

  return (res);

} /* te_load_dct() */


/* Apply the DCT function to a single frame of test data, any DCT routine excepting feature
 * rich fixed point DCTs and blockwise DCTs. */
int te_frameProc_stnd_dct( tTestEngContext * context, 
                     const fract16         * in,
                           float64_t       * out,
                     tVec * xVec, tVec * yVec )
{
  tHandleTbl dct_handle;
  dct_fxnsf_t fxn = NULL;
  tTestEngContext_dct * context_dct;
  void *px, *py;

  int bexp, shift;
  int N, logN;
  int noReuse, doInplace, isRealForward;

  uint32_t crcSum = 0;

  NASSERT( context && context->desc && context->target.handle );
  NASSERT( in && out && xVec && yVec );
  NASSERT( 0 == ( context->args.dim[0] & ( context->args.dim[0] - 1 ) ) );
  NASSERT( 0 == ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) );
  
  context_dct = (tTestEngContext_dct *)context->target.handle;
  N           = context->args.dim[0];
  logN        = 30 - S_exp0_l( N );

  /* If the DCT routine supports inplace operation, try it every second frame. */
  doInplace = ( context->desc->extraParam & TE_FFT_OPT_INPLACE ) && ( context_dct->frameCnt & 1 );
  /* Also check if the target DCT is allowed to reuse the input buffer for intermediate data. */
  noReuse = !( context->desc->extraParam & TE_FFT_OPT_REUSE_INBUF ) && !doInplace;
  /* Real-valued forward DCT requires special block exponent on input. */
  isRealForward = ( context->desc->extraParam & TE_DCT ) &&
                  ( context->args.caseType == TE_FFT_TESTCASE_FORWARD );

  /* For all fixed point DCTs, block exponent of input data should be at least 1, but for
   * real-valued forward DCT zero is allowed. */
  bexp = ( isRealForward ? 0 : 1 );

  /* Convert 16-bit PCM input data to target DCT format. */
  shift = vecFromPcm16( xVec, (fract16*)in, bexp );

  /* Select in/out buffers for the DCT, and wipe the output buffer. */
  if ( doInplace )
  {
    if ( vecGetSize( xVec ) < vecGetSize( yVec ) )
    {
      memcpy( vecGetElem( yVec, 0 ), vecGetElem( xVec, 0 ), vecGetSize( xVec ) );
      px = py = vecGetElem( yVec, 0 );
    }
    else
    {
      px = py = vecGetElem( xVec, 0 );
    }
  }
  else
  {
    memset( vecGetElem( yVec, 0 ), 0, vecGetSize( yVec ) );
    px = vecGetElem( xVec, 0 );
    py = vecGetElem( yVec, 0 );
  }

  /* Select the target DCT routine (either forward or inverse). */
  if ( context->args.caseType == TE_FFT_TESTCASE_FORWARD )
  {
      findItem(&dct_handle,N,(dct_fxns_t)((const tTestEngDesc_dct *)context->desc)->frwTransFxn);
      fxn = (dct_fxnsf_t)dct_handle.fun;
    /* Compensate for scaling shift performed by fixed point DCTs. */
    shift -= logN;
  }
  else if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
  {
      findItem(&dct_handle,N,(dct_fxns_t)((const tTestEngDesc_dct *)context->desc)->invTransFxn);
      fxn = (dct_fxnsf_t)dct_handle.fun;
    /* For fixed point inverse DCT, we have to divide output signal by DCT size.
     * Just don't compensate for the scaling shift performed by the DCT routine. */
  }
  else
  {
    NASSERT( !"Bad test case type!" );
  }

  /* If reuse is not allowed, make sure the buffer stays intact. */
  if ( noReuse ) crcSum = crc32( 0, (uint8_t*)px, vecGetSize( xVec ) );
    /* add to the log */
  {
    char str[30];
    tReportFUT fut[1];
    const tTestEngContext_dct_int * context_dct = (tTestEngContext_dct_int *)context->target.handle;
    fut[0]=(tReportFUT)fxn;
    sprintf(str,"N=%d",dct_handle.N);
    vReportAdd(fut,1,str,context_dct->ext.fInName,context->args.caseType,vecGetSize( xVec )+vecGetSize( yVec ));
  }
  /* Apply the target DCT routine. */
  fxn( py, px, dct_handle.h );

  if ( doInplace && vecGetSize( xVec ) >= vecGetSize( yVec ) )
  {
    memcpy( vecGetElem( yVec, 0 ), py, vecGetSize( yVec ) );
  }
  
  if ( noReuse && crcSum != crc32( 0, (uint8_t*)px, vecGetSize( xVec ) ) )
  {
    printf( "te_frameProc_stnd_dct(): target DCT has corrupted the input buffer\n" );
    return ( 0 );
  }

  /* Convert output data to complex 64-bit floating point and rescale them. */
  vecToFp64( (float64_t*)out, yVec, shift );

  return (1);

} /* te_frameProc_stnd_dct() */

/* Apply the DCT function to a single frame of test data, feature rich fixed point
 * DCT (with scaling method option). */
int te_frameProc_stnd_scl_dct( tTestEngContext * context, 
                         const fract16         * in,
                               float64_t       * out,
                         tVec * xVec, tVec * yVec )
{
  tHandleTbl dct_handle;
  dct_fxns_t fxn = NULL;
  tTestEngContext_dct * context_dct;
  void *px, *py;

  int bexp=0, shift;
  int N, logN;
  int noReuse, doInplace;

  uint32_t crcSum = 0;

  NASSERT( context && context->desc && context->target.handle );
  NASSERT( in && out && xVec && yVec );
  NASSERT( 0 == ( context->args.dim[0] & ( context->args.dim[0] - 1 ) ) );
  NASSERT( 0 != ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) );
  
  context_dct = (tTestEngContext_dct *)context->target.handle;
  N           = context->args.dim[0];
  logN        = 30 - S_exp0_l( N );

  /* If the DCT routine supports inplace operation, try it every second frame. */
  doInplace = ( context->desc->extraParam & TE_FFT_OPT_INPLACE ) && ( context_dct->frameCnt & 1 );
  /* Also check if the target DCT is allowed to reuse the input buffer for intermediate data. */
  noReuse = !( context->desc->extraParam & TE_FFT_OPT_REUSE_INBUF ) && !doInplace;

  /* Select the target DCT routine. */
  if ( context->args.caseType == TE_FFT_TESTCASE_FORWARD )
  {
      findItem(&dct_handle,N,(dct_fxns_t)((const tTestEngDesc_dct *)context->desc)->frwTransFxn);
      fxn = dct_handle.fun;
  }
  else if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
  {
      findItem(&dct_handle,N,(dct_fxns_t)((const tTestEngDesc_dct *)context->desc)->invTransFxn);
      fxn = dct_handle.fun;
  }
  else
  {
    NASSERT( !"Bad test case type!" );
  }

  /* Select the block exponent for fixed point DCT input data. Scale method 0 is
   * "no scaling", 3 - "static scaling". */
  if (context_dct->scale_method == 3 || context_dct->scale_method == 1 || context_dct->scale_method == 2)
  {
      bexp = 0;
  }
  else
  {
      bexp = logN;
  }

  /* Convert 16-bit PCM input data to target DCT format. */
  shift = vecFromPcm16( xVec, (fract16*)in, bexp );

  /* Select in/out buffers for the DCT, and wipe the output buffer. */
  if ( doInplace )
  {
    if ( vecGetSize( xVec ) < vecGetSize( yVec ) )
    {
      memcpy( vecGetElem( yVec, 0 ), vecGetElem( xVec, 0 ), vecGetSize( xVec ) );
      px = py = vecGetElem( yVec, 0 );
    }
    else
    {
      px = py = vecGetElem( xVec, 0 );
    }
  }
  else
  {
    memset( vecGetElem( yVec, 0 ), 0, vecGetSize( yVec ) );
    px = vecGetElem( xVec, 0 );
    py = vecGetElem( yVec, 0 );
  }

  /* If reuse is not allowed, make sure the buffer stays intact. */
  if ( noReuse ) crcSum = crc32( 0, (uint8_t*)px, vecGetSize( xVec ) );

  /* Add to the log. */
  {
    char str[30];
    tReportFUT fut[1];
    //const tTestEngContext_dct_int * context_dct = (tTestEngContext_dct_int *)context->target.handle;
    fut[0]=(tReportFUT)fxn;
    sprintf(str,"N=%d,scale_method=%d",dct_handle.N,context_dct->scale_method);
    vReportAdd(fut,1,str,context_dct->fInName,context->args.caseType,vecGetSize( xVec )+vecGetSize( yVec ));
  }

  /* Apply the target DCT routine. */
  bexp=fxn( py,px,dct_handle.h,context_dct->scale_method); 

  if ( doInplace && vecGetSize( xVec ) >= vecGetSize( yVec ) )
  {
    memcpy( vecGetElem( yVec, 0 ), py, vecGetSize( yVec ) );
  }
  
  if ( noReuse && crcSum != crc32( 0, (uint8_t*)px, vecGetSize( xVec ) ) )
  {
    printf( "te_frameProc_stnd_scl_dct(): target DCT has corrupted the input buffer\n" );
    return ( 0 );
  }

  /* Compensate for scaling shift performed by fixed point forward DCTs. */
  shift -= bexp;
  /* For inverse DCTs, we also have to divide output signal by FFT size. */
  if (context->args.caseType == TE_FFT_TESTCASE_INVERSE)
  {
      shift += logN;
  }

  /* Convert output data to complex 64-bit floating point and rescale them. */
  vecToFp64( (float64_t*)out, yVec, shift );

  return (1);

} /* te_frameProc_stnd_scl_dct() */

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Test data file processing for a standalone test. Applies the target DCT function
 * to test data, compares the DCT output with reference data and estimates the SINAD. */
void fileProc_stnd( tTestEngContext * context )
{
  tTestEngContext_dct           * context_dct;
  tTestEngFrameProcFxn_fft_stnd * frameProcFxn;

  tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
  tVec xVec, yVec;            /* In/out vectors for the target DCT routine. */
  FILE *fIn = 0, *fRef = 0;
  int res = 0;
  int lenIn=0, lenOut=0;
  int fmtX=0, fmtY=0;
  int n, N;

  NASSERT( context && context->desc );
  NASSERT( context->target.fut && context->target.handle );

  context_dct = (tTestEngContext_dct*)context->target.handle;

  /* DCT size */
  N = context->args.dim[0];

  /* Select length and format (real/complex) for all vectors. */
  switch ( context->args.caseType )
  {
  case TE_FFT_TESTCASE_FORWARD:
    lenIn  = N;
    lenOut = N;
    if (context->desc->extraParam & TE_MDCT) lenIn = 2*N;
    break;
  case TE_FFT_TESTCASE_INVERSE:
    lenIn  = N;
    lenOut = N;
    if (context->desc->extraParam & TE_MDCT) lenOut = 2*N;
    break;
  default:
    ASSERT( 0 );
  }
  fmtX   = FMT_REAL;
  fmtY   = FMT_REAL;

  /* Preset an invalid SINAD for the test result. */
  *vecGetElem_fl32( &context->dataSet.Z, 0 ) = -HUGE_VALF;

  if ( context->isVerbose )
  {
    if ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH )
    {
      printf( "scale_mtd %d ", context_dct->scale_method );
    }

    printf( "%-40s  ", context_dct->fRefName );
  }

  memset( &inVec , 0, sizeof(inVec ) );
  memset( &outVec, 0, sizeof(outVec) );
  memset( &refVec, 0, sizeof(refVec) );
    {
        char * fname;
        const char *dir;
        dir=context->isFull ? FULL_VECTOR_DIR :BRIEF_VECTOR_DIR;
        fname=(char*)malloc(strlen(context_dct->fInName)+strlen(dir)+2);
        sprintf(fname,"%s/%s",dir,context_dct->fInName);
        fIn = fopen(fname, "rb");
        free(fname);
        fname=(char*)malloc(strlen(context_dct->fRefName)+strlen(dir)+2);
        sprintf(fname,"%s/%s",dir,context_dct->fRefName);
        fRef = fopen(fname, "rb");
        free(fname);
    }

  /* Allocate vectors for in/out/reference data. */
  if ( !vecAlloc( &inVec , lenIn , TE_ALIGN_YES, FMT_FRACT16|fmtX, 0 ) ||
       !vecAlloc( &outVec, lenOut, TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) ||
       !vecAlloc( &refVec, lenOut, TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) )
  {
    printf( "fileProc_stnd(): failed to allocate in/out/ref vectors, N=%d\n", N );
  }
  /* Open input data file. */
  else if (fIn==NULL )
  {
    printf( "fileProc_stnd(): failed to open %s for reading\n", context_dct->fInName );
  }
  /* Open reference data file. */
  else if (fRef == NULL )
  {
    printf( "fileProc_stnd(): failed to open %s for reading\n", context_dct->fRefName );
  }
  else res = 1;

  /*
   * Process the input file frame-by-frame.
   */

  if ( res )
  {
    int baseFmt = ( context->desc->fmt & FMT_DTYPE_MASK );
    int isAligned = context->desc->isAligned;

    fract16   * pin  = (fract16   *)vecGetElem( &inVec , 0 );
    float64_t * pout = (float64_t *)vecGetElem( &outVec, 0 );
    float64_t * pref = (float64_t *)vecGetElem( &refVec, 0 );

    float32_t sinadAvg, sinadMin = HUGE_VALF;
    float64_t errSum = 0, refSum = 0;

    int efbMin = 32;

    memset( &xVec, 0, sizeof(xVec) );
    memset( &yVec, 0, sizeof(yVec) );

    context_dct->frameCnt = 0;

    frameProcFxn = ( (tTestEngFrameProcFxn_fft_stnd*)context->target.fut );

    /* Read 16-bit complex samples from the input file. */
    while ( ( n = fread( pin, inVec.szElem, ( lenIn & ~1 ), fIn ) ) > 0 )
    {
      /* Augment the (last) incomplete frame with zeros. */
      memset( (uint8_t*)pin + n*inVec.szElem, 0, ( ( lenIn & ~1 ) - n )*inVec.szElem );
      /* Zero the output frame. */
      memset( pout, 0, lenOut*outVec.szElem );

      /* Allocate in/out buffers for the target DCT routine. */
      if ( ( !xVec.szBulk && !vecAlloc( &xVec, lenIn , isAligned, baseFmt | fmtX, 0 ) ) ||
           ( !yVec.szBulk && !vecAlloc( &yVec, lenOut, isAligned, baseFmt | fmtY, 0 ) ) )
      {
        printf( "fileProc_stnd(): failed to allocate xVec/yVec, "
                "frameCnt=%d, N=%d\n", context_dct->frameCnt, N );
        res = 0; break;
      }

      /* Use a proprietary frame processing function to perform the target DCT. */
      if ( !( res = frameProcFxn( context, pin, pout, &xVec, &yVec ) ) ) break;

      /* When in unaligned mode, in/out vectors should be periodically reallocated to try various
       * address offsets. */
      if ( !( ++context_dct->frameCnt % XY_VEC_REALLOC_PERIOD ) && !isAligned )
      {
        if ( !( res = vecsFree( &xVec, &yVec, 0 ) ) )
        {
          printf( "fileProc_stnd(): vecsFree() failed for xVec/yVec, "
                  "frameCnt=%d, N=%d\n", context_dct->frameCnt, N );
          break;
        }

        memset( &xVec, 0, sizeof(xVec) );
        memset( &yVec, 0, sizeof(yVec) );
      }

      /* Read double precision reference samples */
      if ( (int)fread( pref, refVec.szElem, lenOut, fRef ) < lenOut )
      {
        printf( "fileProc_stnd(): failed to read reference data, "
                "frameCnt=%d, N=%d\n", context_dct->frameCnt, lenOut );
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
      sinadMin = 10.f*STDLIB_MATH(log10f)( sinadMin );

      /* Set the test result for test verification. */
      *vecGetElem_fl32( &context->dataSet.Z, 0 ) = sinadMin;

      if ( context->isVerbose )
      {
        sinadAvg = ( refSum > 0 ? 10.f*STDLIB_MATH(log10f)( (float32_t)(refSum/errSum) ) : -HUGE_VALF );

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

/* Allocate vectors and load the data set for dct2d:
 * X[N*N][L], Z[N*N][L] */
int te_loadFxn_dct2d( tTestEngContext * context )
{
    tTestEngContext_dct_int * context_dct = (tTestEngContext_dct_int *)context->target.handle;
    tVec X, Z, Zlo, Zhi;
    int N, L;
    int fmtX, fmtZ;
    int nElem, res = 0;
    fmtX = fmtZ = 0;
    ASSERT( context && context->seqFile );

    memset( &X  , 0, sizeof(X  ) );
    memset( &Z  , 0, sizeof(Z  ) );
    memset( &Zlo, 0, sizeof(Zlo) );
    memset( &Zhi, 0, sizeof(Zhi) );

    N = MAX( 0, context->args.dim[0] );
    L = MAX( 0, context->args.dim[1] );

    nElem = N*N*L;
    if (context->desc->extraParam & TE_DCT2D)
    {
      fmtX = FMT_UINT8|FMT_REAL;
      fmtZ = FMT_INT16|FMT_REAL;
    }
    else if (context->desc->extraParam & TE_IDCT2D)
    {
      fmtX = FMT_INT16|FMT_REAL;
      fmtZ = FMT_UINT8|FMT_REAL;
    }
    else
    {
      ASSERT("Unsupported type of DCT!\n");
    }

    /* If DCT supports the scaling option, read the scaling method from the SEQ-file. */
    if ( ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) &&
         ( 1 != seqFileScanf( context->seqFile, "%d", &context_dct->ext.scale_method ) ) )
    {
      printf( "te_load_dct2d(): bad SEQ-file format\n" );
    }
    else
    {
      res = 1;
    }
    /* Allocate data vectors memory. */
    res&=(0 != vecAlloc( &X, nElem, context->desc->isAligned, fmtX, NULL ));
    res&=(0 != vecAlloc( &Z, nElem, context->desc->isAligned, fmtZ, NULL ));
    res&=(2 == vecsAlloc( TE_ALIGN_NO,
                            fmtZ,
                            &Zlo, nElem,
                            &Zhi, nElem, 0 ));
    if (!res)
    {
    printf( "loadFxn_dct2d(): failed to allocate vectors; "
            "fmtX = 0x%02x, fmtZ = 0x%02x, nElemX = %d, nElemZ = %d\n",
            (unsigned)fmtX, (unsigned)fmtZ, nElem, nElem );
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile, &X, &Zlo, &Zhi, 0 ) )
    {
    printf( "loadFxn_dct2d(): failed to read vectors data; "
            "fmt = 0x%02x, nElemX = %d, nElemZ = %d\n",
            (unsigned)context->desc->fmt, nElem, nElem );
    }
    else
    {
        memset( &context->dataSet, 0, sizeof( context->dataSet ) );

        context->dataSet.X   = X;
        context->dataSet.Z   = Z;
        context->dataSet.Zlo = Zlo;
        context->dataSet.Zhi = Zhi;
        res = 1;
    }

  if ( !res ) te_freeVectors(context); /* Free vectors data if failed. */
  return (res);

} /* loadFxn_dct2d() */

/* Apply the function under test to test case data set: dct2d */
void te_processFxn_dct2d( tTestEngContext * context )
{
  tHandleTbl dct_handle;
  dct_fxns_2d_t fxn = NULL;
  tTestEngContext_dct * context_dct;
  int N, L;
  void *X, *Y;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Z, 0 );

  context_dct = (tTestEngContext_dct *)context->target.handle;
  N = context->args.dim[0];
  L = context->args.dim[1];

  if (context->desc->extraParam & TE_DCT2D)
  {
    findItem(&dct_handle,N,(dct_fxns_t)((const tTestEngDesc_dct *)context->desc)->frwTransFxn);
    fxn = (dct_fxns_2d_t)dct_handle.fun;
  }
  else if (context->desc->extraParam & TE_IDCT2D)
  {
    findItem(&dct_handle,N,(dct_fxns_t)((const tTestEngDesc_dct *)context->desc)->invTransFxn);
    fxn = (dct_fxns_2d_t)dct_handle.fun;
  }
  else
  {
    ASSERT(0);
  }

  te_vReportStd(context);
  /* Apply the target DCT routine. */
  fxn( Y, X, dct_handle.h, L, context_dct->scale_method );

} /* processFxn_dct2d() */
