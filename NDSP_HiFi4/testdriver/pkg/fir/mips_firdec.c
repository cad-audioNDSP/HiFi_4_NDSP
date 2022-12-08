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
* Test module for testing cycle performance (Decimation)
*/

#include "config.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define FIRDEC16X16_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec16x16  ,(D,M),(objinstance_memory, D,M, mips.inp2.i16),(firdec16x16  ,mips.out2.i16, mips.inp2.i16,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC32X32_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec32x32  ,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firdec32x32  ,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC32X32EP_PROFILE(cond,verb,D,N,M) OBJ_PROFILE_INVERTED(cond,verb,firdec32x32ep,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firdec32x32ep,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC32X16_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec32x16  ,(D,M),(objinstance_memory, D,M, mips.inp2.i16),(firdec32x16  ,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC24X24_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec24x24  ,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firdec24x24  ,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDECF_PROFILE(cond,verb,D,N,M)       OBJ_PROFILE_INVERTED(cond,verb,firdecf      ,(D,M),(objinstance_memory, D,M, mips.inp2.f32),(firdecf      ,mips.out2.f32, mips.inp2.f32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);

void mips_firdec16x16(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC16X16_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2,   80, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC16X16_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC16X16_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 7, 1024, 260)
}
void mips_firdec32x16(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 2)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 4 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 8 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 16)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 32)
    FIRDEC32X16_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2,   80, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,  4 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRDEC32X16_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRDEC32X16_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 80, 256)
}
void mips_firdec24x24(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
#if 0// for HiFi3/3z
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC24X24_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC24X24_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC24X24_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2,   80, 256)
#endif
}
void mips_firdec32x32(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,   8)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,  16)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRDEC32X32_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,  4 )
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRDEC32X32_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRDEC32X32_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2,   80, 256)
   FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2,   80, 256)
   FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,   2)
   FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  4 )
   FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  8 )
   FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  16)
   FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRDEC32X32EP_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  4 )
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRDEC32X32EP_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRDEC32X32EP_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 80, 256)
}

#define PROFILE_DECIMA(cond,verb,fun,N,M,D,firstate_type,suffix)                                                              \
{                                                                                                                   \
    firstate_type state;                                                                                            \
    fir_init((&state),mips.inp2. suffix,mips.inp1. suffix,M,D);                                                               \
    PROFILE_INVERTED_NEW(cond,verb,fun,(&state,mips.out0. suffix,mips.inp0. suffix,N),fout,"N=" #N ",M=" #M",D=" #D,prf_maccycle,(M*(N/D)));\
}

#define PROFILE_DECIMAf(cond,verb,fun,N,M,D) PROFILE_DECIMA(cond,verb,fun,N,M,D,fir_statef,f32) 

#define PROFILE_INTERP(cond,verb,fun,N,M,D,firstate_type,suffix)                                                              \
{                                                                                                                   \
    firstate_type state;                                                                                            \
    fir_init((&state),mips.inp2. suffix,mips.inp1. suffix,M,D);                                                               \
    PROFILE_INVERTED_NEW(cond,verb,fun,(&state,mips.out0. suffix,mips.inp0. suffix,N),fout,"N=" #N ",M=" #M",D=" #D,prf_maccycle,(M*N*D));  \
}

#define PROFILE_INTERPf(cond,verb,fun,N,M,D)         PROFILE_INTERP(cond,verb,fun,N,M,D,fir_statef    ,f32)

void mips_firdec(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {
      mips_firdec16x16(phaseNum, isFull, isVerbose, fout);
      mips_firdec32x16(phaseNum, isFull, isVerbose, fout);
      mips_firdec24x24(phaseNum, isFull, isVerbose, fout);
      mips_firdec32x32(phaseNum, isFull, isVerbose, fout);
    }
    if ( phaseNum == 0 || phaseNum == 2 )
    {
        FIRDECF_PROFILE(     1, isVerbose, 2,  1024, 256);
        FIRDECF_PROFILE(isFull, isVerbose, 2,  1024, 512);
        FIRDECF_PROFILE(     1, isVerbose, 3,  1024, 256);
        FIRDECF_PROFILE(isFull, isVerbose, 3,  1024, 512);
        FIRDECF_PROFILE(     1, isVerbose, 4,  1024, 256);
        FIRDECF_PROFILE(isFull, isVerbose, 4,  1024, 512);
        FIRDECF_PROFILE(isFull, isVerbose, 8,  1024, 256);
        FIRDECF_PROFILE(isFull, isVerbose, 8,  1024, 512);
        FIRDECF_PROFILE(isFull, isVerbose, 11, 1024, 256 );
        FIRDECF_PROFILE(isFull, isVerbose, 11, 1024, 512 );
        FIRDECF_PROFILE(isFull, isVerbose, 23, 1024, 256 );
        FIRDECF_PROFILE(isFull, isVerbose, 23, 1024, 512 );
    }
} /* mips_firdec() */
