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
* Test module for testing cycle performance (Scalar Math)
*/

#include "mips.h"
#include "config.h"
#include LIBRARY_HEADER(math)

void mips_maths(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {

        PROFILE_SIMPLE(isFull, isVerbose, scl_recip16x16,  (7002 ),                    fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_recip32x32,  (7002 ),                    fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_recip24x24,  (-1966),                    fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_recip64x64,  (7002 ),                    fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_divide16x16, (-17621, -29508),           fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_divide32x32, (-1154751292, -1933789767), fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_divide24x24, (1154751292, -1933789767),  fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_divide64x32, (-1154751292, -1933789767), fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_divide64x64, (-1154751292, -1933789767), fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_log2_32x32, (496366179),                 fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_logn_32x32, (496366179),                 fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_log10_32x32, (496366179),                fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_log2_24x24, (496366179),                 fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_logn_24x24, (496366179),                 fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_log10_24x24, (496366179),                fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilog2_32x32, (-1010430329),           fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilogn_32x32, (-1010430329),           fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilog10_32x32, (-1010430329),          fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilog2_24x24, (-1010430329),           fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilogn_24x24, (-1010430329),           fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilog10_24x24, (-1010430329),          fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt16x16, (-29508),                     fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt32x16, (-1154751292),                fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt32x32, (-1154751292),                fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt24x24, (-1154751292),                fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt64x32, (-1105678961268363246),       fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_sine32x32, (-1154751292),                fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_cosine32x32, (-1154751292),              fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_sine24x24, (-1154751292),                fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_cosine24x24, (-1154751292),              fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_tan32x32, (2147483640),                  fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_tan24x24, (2147483640),                  fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_atan32x32, (-1154751292),                fout,"",prf_cycle);
#if 0 //HiFi3/3z API
        PROFILE_SIMPLE(isFull, isVerbose, scl_atan24x24, (-1154751292),                fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_atan2_24x24, (-1154751292, -1010430329), fout,"",prf_cycle);
#endif
        PROFILE_SIMPLE(isFull, isVerbose, scl_rsqrt16x16, (29508),                    fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_rsqrt32x32, (1154751292),               fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_sigmoid32x32, (-1154751292),             fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_tanh32x32, (-1154751292),                fout,"",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_relu32x32, (1154751292, 1010430329),     fout,"",prf_cycle);

    }
    if ( phaseNum == 0 || phaseNum == 2 )
    {
        //typedef union {struct { float32_t re, im;}s; complex_float z;} tcomplex_float;

        PROFILE_SIMPLE(isFull, isVerbose, scl_int2float,(9621325  ,1),fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_float2int,(9621325.f,1),fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_sinef     ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_cosinef   ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_tanf      ,(0.4f),      fout,"x=0.4",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_tanf      ,(1.2f),      fout,"x=1.2",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_log2f     ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_log10f    ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_lognf     ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilog2f ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilog10f,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_antilognf ,(1.2f),      fout,""     ,prf_cycle);
        PROFILE_SIMPLE( isFull,isVerbose, scl_powf       ,(1.f, 1.f), fout, "x=1 y=1"        , prf_cycle);
        PROFILE_SIMPLE( isFull,isVerbose, scl_powf       ,(1.25f, 0.75f), fout, "x=1.25 y=0.75"  , prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_atanf     ,(0.7f),      fout,"x=0.7",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_atanf     ,(1.3f),      fout,"x=1.3",prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_atan2f,    (1.2f,2.f),  fout,""     ,prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_sigmoidf,  (1.2f),      fout,"",     prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_tanhf,     (1.2f),      fout,"",     prf_cycle);
        PROFILE_SIMPLE(isFull, isVerbose, scl_reluf,     (1.2f,2.f),  fout,""     ,prf_cycle);
    }
} /* mips_maths() */
