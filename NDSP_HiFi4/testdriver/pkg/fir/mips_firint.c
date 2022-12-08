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
* Test module for testing cycle performance (Interpolation)
*/

#include "config.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define FIRINTERP16X16_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firinterp16x16  ,(D,M),(objinstance_memory, D,M, mips.inp2.i16),(firinterp16x16  ,mips.out2.i16, mips.inp2.i16, N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M*D);
#define FIRINTERP32X32_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firinterp32x32  ,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firinterp32x32  ,mips.out2.i32, mips.inp2.i32, N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M*D);
#define FIRINTERP32X32EP_PROFILE(cond,verb,D,N,M) OBJ_PROFILE_INVERTED(cond,verb,firinterp32x32ep,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firinterp32x32ep,mips.out2.i32, mips.inp2.i32, N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M*D);
#define FIRINTERP32X16_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firinterp32x16  ,(D,M),(objinstance_memory, D,M, mips.inp2.i16),(firinterp32x16  ,mips.out2.i32, mips.inp2.i32, N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M*D);
#define FIRINTERP24X24_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firinterp24x24  ,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firinterp24x24  ,mips.out2.i32, mips.inp2.i32, N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M*D);
#define FIRINTERPF_PROFILE(cond,verb,D,N,M)       OBJ_PROFILE_INVERTED(cond,verb,firinterpf      ,(D,M),(objinstance_memory, D,M, mips.inp2.f32),(firinterpf      ,mips.out2.f32, mips.inp2.f32, N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M*D);

void mips_firinterp16x16(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRINTERP16X16_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 3, 1024,   4)
    FIRINTERP16X16_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRINTERP16X16_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRINTERP16X16_PROFILE(isFull, isVerbose, 2,   80, 204)
}
void mips_firinterp32x16(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 2, 1024,  8 )
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 2, 1024,  16)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRINTERP32X16_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 3, 1024,   4)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 3, 1024,   8)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRINTERP32X16_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRINTERP32X16_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRINTERP32X16_PROFILE(isFull, isVerbose, 8, 80, 204)
}

void mips_firinterp24x24(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
#if 0// for HiFi3/3z
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRINTERP24X24_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 3, 1024,   4)
    FIRINTERP24X24_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRINTERP24X24_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRINTERP24X24_PROFILE(isFull, isVerbose, 2,   80, 204)
#endif
}
void mips_firinterp32x32(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 2, 1024,  8 )
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 2, 1024,  16)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRINTERP32X32_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 3, 1024,   4)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRINTERP32X32_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRINTERP32X32_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRINTERP32X32_PROFILE(isFull, isVerbose, 8, 80, 204)
  
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 2, 80, 204)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  8 )
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  16)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRINTERP32X32EP_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 3, 1024,   4)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRINTERP32X32EP_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRINTERP32X32EP_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRINTERP32X32EP_PROFILE(isFull, isVerbose, 8, 80, 204)
}
#define PROFILE_INTERP(cond,verb,fun,N,M,D,firstate_type,suffix)                                                              \
{                                                                                                                   \
    firstate_type state;                                                                                            \
    fir_init((&state),mips.inp2. suffix,mips.inp1. suffix,M,D);                                                               \
    PROFILE_INVERTED_NEW(cond,verb,fun,(&state,mips.out0. suffix,mips.inp0. suffix,N),fout,"N=" #N ",M=" #M",D=" #D,prf_maccycle,(M*N*D));  \
}

#define PROFILE_INTERPf(cond,verb,fun,N,M,D)         PROFILE_INTERP(cond,verb,fun,N,M,D,fir_statef    ,f32)

void mips_firint(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {
      mips_firinterp16x16(phaseNum, isFull, isVerbose, fout);
      mips_firinterp32x16(phaseNum, isFull, isVerbose, fout);
      mips_firinterp24x24(phaseNum, isFull, isVerbose, fout);
      mips_firinterp32x32(phaseNum, isFull, isVerbose, fout);
    }
    if ( phaseNum == 0 || phaseNum == 2 )
    {
        FIRINTERPF_PROFILE(     1, isVerbose, 2,1024,256);
        FIRINTERPF_PROFILE(isFull, isVerbose, 2,1024,512);
        FIRINTERPF_PROFILE(     1, isVerbose, 3,1024,256);
        FIRINTERPF_PROFILE(isFull, isVerbose, 3,1024,512);
        FIRINTERPF_PROFILE(     1, isVerbose, 4,1024,256);
        FIRINTERPF_PROFILE(isFull, isVerbose, 4,1024,512);
        FIRINTERPF_PROFILE(isFull, isVerbose, 8,1024,256);
        FIRINTERPF_PROFILE(isFull, isVerbose, 8,1024,512);
    }
} /* mips_firint() */
