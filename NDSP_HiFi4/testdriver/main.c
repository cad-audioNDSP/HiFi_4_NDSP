/* ------------------------------------------------------------------------ */
/* Copyright (c) 2018 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ("Cadence    */
/* Libraries") are the copyrighted works of Cadence Design Systems Inc.     */
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
  Test module/example for testing NatureDSP_Signal library
  Integrit, 2006-2018
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/* Portable data types */
#include "types.h"
/* Fixed-point arithmetic primitives. */
#include "NatureDSP_Math.h"
/* Test environmant configuration */
#include "config.h"
/* Cycles measurement API. */
#include "mips.h"
#include "vectools.h"
#include "vreport.h"

int main_firblk   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_firblk   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_firdec   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_firdec   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_firint   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_firint   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_firother ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_firother ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_conv2d   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_conv2d   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_iirbq    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_iirbq    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_iirbq_nd    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_iirbq_nd    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_iirlt    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_iirlt    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_mathv    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_mathv    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_mathvf   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_mathvf   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_maths    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_maths    ( int phaseNum, int isFull, int isVerbose, FILE* f ); int accuracy_maths(int phaseNum, int isFull, int isVerbose, int breakOnError, int optAccuracy, int optException);
int main_complexv ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_complexv ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_complexs ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_complexs ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_vector   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_vector   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_ef       ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_ef       ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_matop    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_matop    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_matinv   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_matinv   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_chol     ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_chol     ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_pfit     ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_pfit     ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_fft      ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_fft      ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_cfft     ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_cfft     ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_rfft     ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_rfft     ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_cnfft    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_cnfft    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_rnfft    ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_rnfft    ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_cfftie   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_cfftie   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_rfftie   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_rfftie   ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_dct      ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_dct      ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_spectrum ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_spectrum ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_mfcc     ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_mfcc     ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_imgrotate( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_imgrotate( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_imgresize( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_imgresize( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_imgmisc  ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_imgmisc  ( int phaseNum, int isFull, int isVerbose, FILE* f );
int main_imgfft   ( int phaseNum, int isFull, int isVerbose, int breakOnError ); void mips_imgfft   ( int phaseNum, int isFull, int isVerbose, FILE* f );

int main_accuracy (int phaseNum, int isFull, int isVerbose, int breakOnError) { return 0;}
int main_exception(int phaseNum, int isFull, int isVerbose, int breakOnError) { return 0; }
int main_mips     ( int phaseNum, int isFull, int isVerbose, int breakOnError );
int main_func(int phaseNum, int isFull, int isVerbose, int breakOnError) { return 0;}
int main_help     ( int phaseNum, int isFull, int isVerbose, int breakOnError );


typedef struct
{
  const char * cmd;
  const char * help;
  int (*fntest)(int,int,int,int);
  void (*fnmips)(int,int,int,FILE*);
  int (*fnaccuracy)(int,int,int,int,int,int);
  uint64_t mask;
}
tParamTbl;

#define PACKAGE_MASK   0x0003ffffffffffffULL
#define PHASE_MASK    0x0700000000000000ULL
#define MIPS_FLAG     0x0800000000000000ULL
#define NOABORT_FLAG  0x1000000000000000ULL
#define VERBOSE_FLAG  0x2000000000000000ULL
#define VERIFY_FLAG   0x0080000000000000ULL
#define FUNC_FLAG     0x0040000000000000ULL
#define FULL_FLAG     0x4000000000000000ULL
#define BRIEF_FLAG    0x0020000000000000ULL
#define NOWARMUP_FLAG 0x0010000000000000ULL
#define ACCURACY_FLAG  0x0008000000000000ULL
#define EXCEPTION_FLAG 0x0004000000000000ULL
#define HELP_FLAG     0x8000000000000000ULL

static const tParamTbl tbl[]=
{
  { "-fir"          , "FIR Filters"                                    , NULL          , NULL              ,NULL          , 0x000000000000001FULL },
  { "-firblk"       , "Filtering"                                      , main_firblk   , mips_firblk       ,NULL          , 0x0000000000000001ULL },
  { "-firdec"       , "Decimation"                                     , main_firdec   , mips_firdec       ,NULL          , 0x0000000000000002ULL },
  { "-firint"       , "Interpolation"                                  , main_firint   , mips_firint       ,NULL          , 0x0000000000000004ULL },
  { "-firother"     , "Correlation, Convolution, Despreading, LMS"     , main_firother , mips_firother     ,NULL          , 0x0000000000000008ULL },
  { "-conv2d"       , "2D convolution"                                 , main_conv2d   , mips_conv2d       ,NULL          , 0x0000000000000010ULL },
  { "-iir"          , "IIR Filters"                                    , NULL          , NULL              ,NULL          , 0x0000000000000070ULL<<1 },
  { "-iirbq"        , "Biquad Filters"                                 , main_iirbq    , mips_iirbq        ,NULL          , 0x0000000000000010ULL<<1 },
  { "-iirbqnd"      , "Biquad Filters, no delay"                       , main_iirbq_nd , mips_iirbq_nd     ,NULL          , 0x0000000000000020ULL<<1 },//do not forget to change constant ot an apropriate one
  { "-iirlt"        , "Lattice Filters"                                , main_iirlt    , mips_iirlt        ,NULL          , 0x0000000000000020ULL<<2 },
  { "-math"         , "Math Functions"                                 , NULL          , NULL              ,NULL          , 0x00000000000001C0ULL<<2 },
  { "-mathv"        , "Vectorized Math"                                , main_mathv    , mips_mathv        ,NULL          , 0x0000000000000040ULL<<2 }, 
  { "-mathvf"       , "Vectorized Fast Math"                           , main_mathvf   , mips_mathvf       ,NULL          , 0x0000000000000080ULL<<2 }, 
  { "-maths"        , "Scalar Math"                                    , main_maths    , mips_maths        ,accuracy_maths, 0x0000000000000100ULL<<2 }, 
  { "-complex"      , "Complex Functions"                              , NULL          , NULL              ,NULL          , 0x0000000000000600ULL<<2 }, 
  { "-complexv"     , "Vectorized Complex Math"                        , main_complexv , mips_complexv     ,NULL          , 0x0000000000000200ULL<<2 }, 
  { "-complexs"     , "Scalar Complex Math"                            , main_complexs , mips_complexs     ,NULL          , 0x0000000000000400ULL<<2 }, 
  { "-vector"       , "\nVector Operations"                            , main_vector   , mips_vector       ,NULL          , 0x0000000000000800ULL<<2 }, 
  { "-ef"           , "\nEmulated Floating Point Operations"           , main_ef       , mips_ef           ,NULL          , 0x0000000000001000ULL<<2 }, 
  { "-matop"        , "\nMatrix Operations"                            , main_matop    , mips_matop        ,NULL          , 0x0000000000001000ULL<<3 }, 
  { "-matinv"       , "Matrix Decomposition and Inversion"             , NULL          , NULL              ,NULL          , 0x0000000000006000ULL<<3 }, 
  { "-gj"           , "Gauss-Jordan"                                   , main_matinv   , mips_matinv       ,NULL          , 0x0000000000002000ULL<<3 }, 
  { "-chol"         , "Cholesky"                                       , main_chol     , mips_chol         ,NULL          , 0x0000000000004000ULL<<3 }, 
  { "-fit"          , "Fitting and Interpolation"                      , NULL          , NULL              ,NULL          , 0x0000000000040000ULL }, 
  { "-pfit"         , "Polynomial Fitting"                             , main_pfit     , mips_pfit         ,NULL          , 0x0000000000040000ULL }, 
  { "-fft"          , "FFT Routines"                                   , NULL          , NULL              ,NULL          , 0x0000000003F80000ULL },
  { "-cfft"         , "Complex FFT"                                    , main_cfft     , mips_cfft         ,NULL          , 0x0000000000080000ULL },
  { "-rfft"         , "Real FFT"                                       , main_rfft     , mips_rfft         ,NULL          , 0x0000000000100000ULL },
  { "-cnfft"        , "Mixed Radix Complex FFT"                        , main_cnfft    , mips_cnfft        ,NULL          , 0x0000000000200000ULL },
  { "-rnfft"        , "Mixed Radix Real FFT"                           , main_rnfft    , mips_rnfft        ,NULL          , 0x0000000000400000ULL },
  { "-cfftie"       , "Complex FFT with Optimized Memory"              , main_cfftie   , mips_cfftie       ,NULL          , 0x0000000000800000ULL },
  { "-rfftie"       , "Real FFT with Optimized Memory"                 , main_rfftie   , mips_rfftie       ,NULL          , 0x0000000001000000ULL },
  { "-dct"          , "DCT"                                            , main_dct      , mips_dct          ,NULL          , 0x0000000002000000ULL },
  { "-spectrum"     , "FFT power spectrum functions"                   , main_spectrum, mips_spectrum      ,NULL          , 0x0000000004000000ULL },
  { "-mfcc"         , "MFCC features extraction"                       , main_mfcc     , mips_mfcc         ,NULL          , 0x0000000008000000ULL },
  { "-img"          , "image processing functions"                     , NULL          , NULL              ,NULL          , 0x00000000f0000000ULL },
  { "-imgrotate"    , "image rotation"                                 , main_imgrotate, mips_imgrotate    ,NULL          , 0x0000000010000000ULL },
  { "-imgresize"    , "image resize"                                   , main_imgresize, mips_imgresize    ,NULL          , 0x0000000020000000ULL },
  { "-imgmisc"      , "miscellaneous image processing"                 , main_imgmisc  , mips_imgmisc      ,NULL          , 0x0000000040000000ULL },
  { "-imgfft"       , "2D FFT for image data"                          , main_imgfft   , mips_imgfft       ,NULL          , 0x0000000080000000ULL },
  { "-mips"         , "test for performance"                           , main_mips     , NULL              ,NULL          , MIPS_FLAG    },
  { "-phase<n>"     , "limit test coverage to library phase n"         , NULL          , NULL              ,NULL          , PHASE_MASK   },
  { "-noabort"      , "do not stop after the first failure"            , NULL          , NULL              ,NULL          , NOABORT_FLAG },
  { "-verbose"      , "print test progess info"                        , NULL          , NULL              ,NULL          , VERBOSE_FLAG },
  { "-vreport"      , "print verification test report"                 , NULL          , NULL              ,NULL          , VERIFY_FLAG  },
  { "-func"         , "functional testing"                             , main_func     , NULL              ,NULL          , FUNC_FLAG    },
  { "-full"         , "use full-size test vectors"                     , NULL          , NULL              ,NULL          , FULL_FLAG    },
  { "-brief"        , "use shorter test vectors"                       , NULL          , NULL              ,NULL          , BRIEF_FLAG   },
  { "-nowarmup"     , "invalidate cache before each cycle measurement" , NULL          , NULL              ,NULL          , NOWARMUP_FLAG},
  { "-accuracy"     , "testing for accuracy"                           , main_accuracy , NULL              ,NULL          , ACCURACY_FLAG},
  { "-exception"    , "testing for exceptions"                         , main_exception, NULL			   ,NULL          , EXCEPTION_FLAG},
  { "-help"         , "this help"                                      , main_help     , NULL              ,NULL          , HELP_FLAG    },
  { "-h"            , 0                                                , main_help     , NULL              ,NULL          , HELP_FLAG    }
};

int main_help( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int k;

  printf("Available switches:\n");
  for (k=0; k<(int)(sizeof(tbl)/sizeof(tbl[0])); k++) 
  {
      if (tbl[k].help==NULL) continue;
      printf("%-8s  %s\n",tbl[k].cmd,tbl[k].help);
  } 

  return (1);

} /* main_help() */
int main( int argc, char * argv[] )
{
  const int tblSize = sizeof(tbl)/sizeof(tbl[0]);
  uint64_t flags = 0, packageMask;
  int phaseNum;
  int m, n;
  /*----------------------------------------------------------------------------*
   *                            Show Library info                               *
   *----------------------------------------------------------------------------*/

  {
    char lib_version[32], api_version[32];

    GET_LIBRARY_VERSION    ( lib_version );
    GET_LIBRARY_API_VERSION( api_version );

    printf( "%s library version: %s\n"
            "%s API version:     %s\n",
            LIBRARY_NAME, lib_version, LIBRARY_NAME, api_version );
  }

  /*----------------------------------------------------------------------------*
   *                         Scan the command line                              *
   *----------------------------------------------------------------------------*/

  for ( m=1; m<argc; m++ )
  {
    for ( n=0; n<tblSize; n++ )
    {
      size_t len = strcspn( tbl[n].cmd, "<[" );

      if ( len < strlen( tbl[n].cmd ) )
      {
        /* Command line option with numeric argument. */
        if ( !strncmp( argv[m], tbl[n].cmd, len ) && isdigit( (int)argv[m][len] ) )
        {
          uint64_t lsb, p, mask = tbl[n].mask;
          /* Get 1 in the lowest bit position of a mask, and zeros elsewhere. */
          lsb = ( (mask<<1) ^ mask ) & mask;
          /* Retrieve the actual parameter. */
          p = lsb*( argv[m][len] - '0' );
          /* Check if actual parameter is non-zero and matches the mask.  */
          if ( p>0 && ( p & mask ) == p ) { flags |= p; break; }
        }
      }
      /* Simple command line option with no argument. */
      else if ( !strcmp( argv[m], tbl[n].cmd ) )
      {
        flags |= tbl[n].mask; break;
      }
    }

    if ( n >= tblSize )
    {
      printf( "Invalid command line option: %s\n", argv[m] );
      return (-1);
    }
  }

  /* Retrieve a library phase number, if 0 then test everything */
  {
    /* Get 1 in the lowest bit position of the mask, and zeros elsewhere. */
    const uint64_t lsb = ( (PHASE_MASK<<1) ^ PHASE_MASK ) & PHASE_MASK;
    /* Place the phase nubmer at lowest bit positions. */
    phaseNum = (int)( ( flags & PHASE_MASK ) >> ( 32+30 - S_exp0_l( (uint32_t)(lsb>>32) ) ) );
  }

  noWarmup=(flags & NOWARMUP_FLAG)?1:0;
#if defined (COMPILER_XTENSA) && (MEM_MODEL!=2)
  if (noWarmup)
  {
      printf("-nowarmup option is ignored. It should be used if code is built with MEM_MODEL=2\n");
  }
#endif

  /* If no particular packages are specified, then we test all the packages. */
  packageMask = ( !( flags & PACKAGE_MASK ) ? PACKAGE_MASK : ( flags & PACKAGE_MASK ) );
  if((flags & (MIPS_FLAG|FUNC_FLAG|ACCURACY_FLAG|EXCEPTION_FLAG))==0) flags |= MIPS_FLAG;    /* select testing the cycles if no -func, -mips specified */
  /*----------------------------------------------------------------------------*
   *                         Run the selected action                            *
   *----------------------------------------------------------------------------*/

  if ( flags & HELP_FLAG )
  {
    main_help( 0, 0, 0, 0 );
    return 0;
  }
  if (flags & MIPS_FLAG)
  {
    uint64_t printedFamily=0;
    FILE * fout;
    int optFull    = ( 0 != ( flags & FULL_FLAG    ) );
    int optVerbose = ( 0 != ( flags & VERBOSE_FLAG ) );

    fout = fopen( "performance.txt", "wb" );
    perf_info( fout, "Function Name\tCycles Measurements\t\n\tInvocation Parameters\tCycles\n" );

    /* Initialize mips testing. */
    main_mips( 0, 0, 0, 0 );

    for ( n=0; n<tblSize; n++ )
    {
        if ( 0 != ( tbl[n].mask & packageMask ) && NULL != tbl[n].fnmips )
        {
            int k;
            /* find family name */
            for ( k=0; k<tblSize; k++ )
            {
                uint64_t t;
                t = tbl[k].mask & tbl[n].mask;
                if (t!=0 && tbl[k].fntest==NULL) break;
            }
            if (k!=tblSize)
            {
                    /* print family name first and test name next */
                if ((printedFamily & tbl[k].mask)==0)
                {
                    perf_info(fout, "\n%s\n",tbl[k].help);
                    printedFamily |= tbl[k].mask;   /* do not print family name for all next members */
                }
            }
            /* print test name */
            perf_info(fout,"%s\n",tbl[n].help);
            tbl[n].fnmips( phaseNum, optFull, optVerbose, fout );
        }
    }

    fclose( fout );

    printf("================= test completed =================\n");
  }
  if (flags & FUNC_FLAG)
  {
    uint64_t printedFamily=0;
    int optFull    = ( 0 != ( flags & FULL_FLAG    ) );
    int optVerbose = ( 0 != ( flags & VERBOSE_FLAG ) );
    int optBreak   = ( 0 == ( flags & NOABORT_FLAG ) );
    vReportInit( (flags & VERIFY_FLAG) ? 1:0,argv[0] );

    for ( n=0; n<tblSize; n++ )
    {
        int k;
        if ( 0 == ( tbl[n].mask & packageMask ) || NULL == tbl[n].fntest ) continue;
            for ( k=0; k<tblSize; k++ )
            {
                uint64_t t;
                t = tbl[k].mask & tbl[n].mask;
                if (t!=0 && tbl[k].fntest==NULL) break;
            }
            if (k!=tblSize)
            {
                    /* print family name first and test name next */
                if ((printedFamily & tbl[k].mask)==0)
                {
                    printf("\n%s\n",tbl[k].help);
                    printedFamily |= tbl[k].mask;   /* do not print family name for all next members */
                }
            }
      /* print test name */
      printf("%s\n",tbl[n].help);
      /* enable allocation of temorary data in memory region taken by performance measurement utilities */
      vecInitRegion(ALLOC_DRAM0);
      /* start tests */
      if ( !tbl[n].fntest( phaseNum, optFull, optVerbose, optBreak ) && optBreak ) break;
    }

    printf("================= test completed =================\n");
    vReportPrint( );
    vReportClose( );

  }
  // accuracy tests
  if (flags & (ACCURACY_FLAG|EXCEPTION_FLAG))
  {
    uint64_t printedFamily=0;
    int optFull    = ( 0 != ( flags & FULL_FLAG    ) );
    int optVerbose = ( 0 != ( flags & VERBOSE_FLAG ) );
    int optBreak   = ( 0 == ( flags & NOABORT_FLAG ) );

    for ( n=0; n<tblSize; n++ )
    {
        int k;
        if ( 0 == ( tbl[n].mask & packageMask ) || NULL == tbl[n].fnaccuracy ) continue;
            for ( k=0; k<tblSize; k++ )
            {
                uint64_t t;
                t = tbl[k].mask & tbl[n].mask;
                if (t!=0 && tbl[k].fnaccuracy==NULL) break;
            }
            if (k!=tblSize)
            {
                    /* print family name first and test name next */
                if ((printedFamily & tbl[k].mask)==0)
                {
                    printf("\n%s\n",tbl[k].help);
                    printedFamily |= tbl[k].mask;   /* do not print family name for all next members */
                }
            }
      /* print test name */
      printf("%s\n",tbl[n].help);
	  if (!tbl[n].fnaccuracy(phaseNum, optFull, optVerbose, optBreak,
		  (flags & ACCURACY_FLAG) ? 1 : 0, (flags & EXCEPTION_FLAG) ? 1 : 0)) break;
    }
    printf("================= test completed =================\n");
  }

  return (0);

} /* main() */
