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
* Test module for testing cycle performance (Matrix Operations)
*/
#include "config.h"
#include LIBRARY_HEADER(matop)
#include "mips.h"

#define PROFILE_mpy8x8(cond,verb,M,N,P)        PROFILE_INVERTED(cond,verb,mtx_mpy8x8,     (pScr,mips.out0.i8, mips.inp0.i8,mips.inp1.i8, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy8x8_fast(cond,verb,M,N,P)   PROFILE_INVERTED(cond,verb,mtx_mpy8x8_fast,(pScr,mips.out0.i8, mips.inp0.i8,mips.inp1.i8, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy8x16(cond,verb,M,N,P)       PROFILE_INVERTED(cond,verb,mtx_mpy8x16,     (pScr,mips.out0.i16, mips.inp0.i8,mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy8x16_fast(cond,verb,M,N,P)  PROFILE_INVERTED(cond,verb,mtx_mpy8x16_fast,(pScr,mips.out0.i16, mips.inp0.i8,mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy16x16(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpy16x16,(pScr,mips.out0.i16, mips.inp0.i16,mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy16x16_fast(cond,verb,M,N,P) PROFILE_INVERTED(cond,verb,mtx_mpy16x16_fast,(pScr,mips.out0.i16, mips.inp0.i16, mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy24x24(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpy24x24,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy24x24_fast(cond,verb,M,N,P) PROFILE_INVERTED(cond,verb,mtx_mpy24x24_fast,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy32x32(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpy32x32,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy32x32_fast(cond,verb,M,N,P) PROFILE_INVERTED(cond,verb,mtx_mpy32x32_fast,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyf(cond,verb,M,N,P)          PROFILE_INVERTED(cond,verb,mtx_mpyf,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyf_fast(cond,verb,M,N,P)     PROFILE_INVERTED(cond,verb,mtx_mpyf_fast,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));

#define PROFILE_mpyt8x8(cond,verb,M,N,P)           PROFILE_INVERTED(cond,verb,mtx_mpyt8x8,(pScr,mips.out0.i8, mips.inp0.i8,mips.inp1.i8, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt8x8_fast(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpyt8x8_fast,(pScr,mips.out0.i8, mips.inp0.i8, mips.inp1.i8, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt8x16(cond,verb,M,N,P)          PROFILE_INVERTED(cond,verb,mtx_mpyt8x16,     (pScr,mips.out0.i16, mips.inp0.i8, mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt8x16_fast(cond,verb,M,N,P)     PROFILE_INVERTED(cond,verb,mtx_mpyt8x16_fast,(pScr,mips.out0.i16, mips.inp0.i8, mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt16x16(cond,verb,M,N,P)         PROFILE_INVERTED(cond,verb,mtx_mpyt16x16,(pScr,mips.out0.i16, mips.inp0.i16,mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt16x16_fast(cond,verb,M,N,P)    PROFILE_INVERTED(cond,verb,mtx_mpyt16x16_fast,(pScr,mips.out0.i16, mips.inp0.i16, mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt32x32(cond,verb,M,N,P)         PROFILE_INVERTED(cond,verb,mtx_mpyt32x32,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyt32x32_fast(cond,verb,M,N,P)    PROFILE_INVERTED(cond,verb,mtx_mpyt32x32_fast,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpytf(cond,verb,M,N,P)             PROFILE_INVERTED(cond,verb,mtx_mpytf,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpytf_fast(cond,verb,M,N,P)        PROFILE_INVERTED(cond,verb,mtx_mpytf_fast,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mtx_vecmpy8x8(cond,verb,fun,M,N)   PROFILE_INVERTED(cond,verb,fun,(mips.out0.i8, mips.inp0.i8, mips.inp1.i8, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpy8x16(cond,verb,fun,M,N)  PROFILE_INVERTED(cond,verb,fun,(mips.out0.i16, mips.inp0.i8, mips.inp1.i16, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpy16x16(cond,verb,fun,M,N) PROFILE_INVERTED(cond,verb,fun,(mips.out0.i16, mips.inp0.i16, mips.inp1.i16, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpy24x24(cond,verb,fun,M,N) PROFILE_INVERTED(cond,verb,fun,(mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpy32x32(cond,verb,fun,M,N) PROFILE_INVERTED(cond,verb,fun,(mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpyf(cond,verb,fun,M,N)     PROFILE_INVERTED(cond,verb,fun,(mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_TRANSP(isFull,isVerbose,fun,coef,M,N,suffix)           \
{                                                                      \
    ASSERT(M*N*sizeof(mips.inp0.suffix[0])<=sizeof(mips.inp0));                  \
    ASSERT(M*N*sizeof(mips.out0.suffix[0])<=sizeof(mips.out0));                  \
    PROFILE_INVERTED(isFull,isVerbose,fun,(mips.out0.suffix,mips.inp0.suffix,M,N),fout,"M=" #M ",N=" #N ,prf_ptscycle2,(1*M*N)); \
}

#define MAX_M 40
#define MAX_N 100
#define MAX_P 8

static void mips_matop_mpy(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.i32;
  if (phaseNum == 0 || phaseNum == 1)
  {
    PROFILE_mpy8x8(     1, isVerbose, 16, 16, 16);
    PROFILE_mpy8x8(     1, isVerbose, 32, 32, 32);
    PROFILE_mpy8x8(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpy8x8(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpy8x8(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpy8x8(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpy8x8(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpy8x8(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpy8x8(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpy8x8(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpy8x8(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpy8x8_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpy8x8_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpy8x8_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpy8x8_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpy8x8_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpy8x8_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpy8x8_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpyt8x8(     1, isVerbose, 16, 16, 16);
    PROFILE_mpyt8x8(     1, isVerbose, 32, 32, 32);
    PROFILE_mpyt8x8(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpyt8x8(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpyt8x8(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpyt8x8(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpyt8x8(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpyt8x8(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpyt8x8(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpyt8x8(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpyt8x8(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpyt8x8_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpyt8x8_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpyt8x8_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpyt8x8_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpyt8x8_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpyt8x8_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpyt8x8_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpy8x16(     1, isVerbose, 16, 16, 16);
    PROFILE_mpy8x16(     1, isVerbose, 32, 32, 32);
    PROFILE_mpy8x16(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpy8x16(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpy8x16(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpy8x16(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpy8x16(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpy8x16(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpy8x16(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpy8x16(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpy8x16(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpy8x16_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpy8x16_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpy8x16_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpy8x16_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpy8x16_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpy8x16_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpy8x16_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpyt8x16(     1, isVerbose, 16, 16, 16);
    PROFILE_mpyt8x16(     1, isVerbose, 32, 32, 32);
    PROFILE_mpyt8x16(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpyt8x16(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpyt8x16(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpyt8x16(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpyt8x16(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpyt8x16(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpyt8x16(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpyt8x16(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpyt8x16(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpyt8x16_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpyt8x16_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpyt8x16_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpyt8x16_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpyt8x16_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpyt8x16_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpyt8x16_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpy16x16(     1, isVerbose, 16, 16, 16);
    PROFILE_mpy16x16(     1, isVerbose, 32, 32, 32);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 84, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 85, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 86, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 87, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 40, 88, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpy16x16(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpy16x16(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpy16x16(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpy16x16(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpy16x16_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpy16x16_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpy16x16_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpy16x16_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpy16x16_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpy16x16_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpy16x16_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpyt16x16(1, isVerbose, 16, 16, 16);
    PROFILE_mpyt16x16(1, isVerbose, 32, 32, 32);
    PROFILE_mpyt16x16(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpyt16x16(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpyt16x16(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpyt16x16(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpyt16x16(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpyt16x16(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpyt16x16(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpyt16x16(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpyt16x16(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpyt16x16_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpyt16x16_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpyt16x16_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpyt16x16_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpyt16x16_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpyt16x16_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpyt16x16_fast(isFull, isVerbose, 8, 32, 4);
#if 0// HiFi3/3z API
    PROFILE_mpy24x24(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpy24x24(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpy24x24(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpy24x24(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpy24x24(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpy24x24(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpy24x24(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpy24x24(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpy24x24(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpy24x24_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpy24x24_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpy24x24_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpy24x24_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpy24x24_fast(isFull, isVerbose, 8, 32, 4);
#endif
    PROFILE_mpy32x32(1, isVerbose, 16, 16, 16);
    PROFILE_mpy32x32(1, isVerbose, 32, 32, 32);
    PROFILE_mpy32x32(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpy32x32(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpy32x32(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpy32x32(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpy32x32(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpy32x32(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpy32x32(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpy32x32(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpy32x32(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpy32x32_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpy32x32_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpy32x32_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpy32x32_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpy32x32_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpy32x32_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpy32x32_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpyt32x32(1, isVerbose, 16, 16, 16);
    PROFILE_mpyt32x32(1, isVerbose, 32, 32, 32);
    PROFILE_mpyt32x32(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpyt32x32(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpyt32x32(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpyt32x32(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpyt32x32(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpyt32x32(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpyt32x32(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpyt32x32(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpyt32x32(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpyt32x32_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpyt32x32_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpyt32x32_fast(1, isVerbose, 8, 80, 4);
    PROFILE_mpyt32x32_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpyt32x32_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpyt32x32_fast(isFull, isVerbose, 8, 16, 4);
    PROFILE_mpyt32x32_fast(isFull, isVerbose, 8, 32, 4);


  }
}

static void mips_matop_vecmpy(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (phaseNum == 0 || phaseNum == 1)
  {
    PROFILE_mtx_vecmpy8x8(isFull, isVerbose, mtx_vecmpy8x8, 16, 100);
    PROFILE_mtx_vecmpy8x8(isFull, isVerbose, mtx_vecmpy8x8, 16, 104);
    PROFILE_mtx_vecmpy8x8(isFull, isVerbose, mtx_vecmpy8x8, 40, 40);
    PROFILE_mtx_vecmpy8x8(isFull, isVerbose, mtx_vecmpy8x8_fast, 16, 100);
    PROFILE_mtx_vecmpy8x8(      1, isVerbose, mtx_vecmpy8x8_fast, 16, 104);
    PROFILE_mtx_vecmpy8x8(isFull, isVerbose, mtx_vecmpy8x8_fast, 40, 40);

    PROFILE_mtx_vecmpy8x16(isFull, isVerbose, mtx_vecmpy8x16, 16, 100);
    PROFILE_mtx_vecmpy8x16(isFull, isVerbose, mtx_vecmpy8x16, 16, 104);
    PROFILE_mtx_vecmpy8x16(isFull, isVerbose, mtx_vecmpy8x16, 40, 40);
    PROFILE_mtx_vecmpy8x16(isFull, isVerbose, mtx_vecmpy8x16_fast, 16, 100);
    PROFILE_mtx_vecmpy8x16(     1, isVerbose, mtx_vecmpy8x16_fast, 16, 104);
    PROFILE_mtx_vecmpy8x16(isFull, isVerbose, mtx_vecmpy8x16_fast, 40, 40);

    PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16, 16, 100);
    PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16, 16, 104);
    PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16, 40, 40);
    PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16_fast, 16, 100);
    PROFILE_mtx_vecmpy16x16(      1, isVerbose, mtx_vecmpy16x16_fast, 16, 104);
    PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16_fast, 40, 40);
#if 0// HiFi3/3z API
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24, 16, 100);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24, 16, 101);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24, 16, 102);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24, 16, 103);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24, 16, 104);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24, 40, 40);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24_fast, 16, 100);
    PROFILE_mtx_vecmpy24x24(1, isVerbose, mtx_vecmpy24x24_fast, 16, 104);
    PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24_fast, 40, 40);
#endif
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32, 16, 100);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32, 16, 101);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32, 16, 102);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32, 16, 103);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32, 16, 104);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32, 40, 40);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32_fast, 16, 100);
    PROFILE_mtx_vecmpy32x32(1, isVerbose, mtx_vecmpy32x32_fast, 16, 104);
    PROFILE_mtx_vecmpy32x32(isFull, isVerbose, mtx_vecmpy32x32_fast, 40, 40);

  }
}

static void mips_matop_transp(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (phaseNum == 0 || phaseNum == 1)
  {

    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8, sizeof(int8_t), 16, 16, i8);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8, sizeof(int8_t), 27, 27, i8);
    PROFILE_TRANSP(     1, isVerbose, mtx_transpose8x8, sizeof(int8_t), 32, 32, i8);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8, sizeof(int8_t), 39, 39, i8);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8, sizeof(int8_t), 48, 48, i8);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8_fast, sizeof(int8_t), 8, 8, i8);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8_fast, sizeof(int8_t), 16, 16, i8);
    PROFILE_TRANSP(     1, isVerbose, mtx_transpose8x8_fast, sizeof(int8_t), 32, 32, i8);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose8x8_fast, sizeof(int8_t), 48, 48, i8);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16, sizeof(int16_t), 16, 16, i16);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16, sizeof(int16_t), 27, 27, i16);
    PROFILE_TRANSP(     1, isVerbose, mtx_transpose16x16, sizeof(int16_t), 32, 32, i16);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16, sizeof(int16_t), 39, 39, i16);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16, sizeof(int16_t), 48, 48, i16);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16_fast, sizeof(int16_t), 8, 8, i16);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16_fast, sizeof(int16_t), 16, 16, i16);
    PROFILE_TRANSP(     1, isVerbose, mtx_transpose16x16_fast, sizeof(int16_t), 32, 32, i16);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose16x16_fast, sizeof(int16_t), 48, 48, i16);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32, sizeof(int32_t), 16, 16, i32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32, sizeof(int32_t), 27, 27, i32);
    PROFILE_TRANSP(     1, isVerbose, mtx_transpose32x32, sizeof(int32_t), 32, 32, i32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32, sizeof(int32_t), 39, 39, i32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32, sizeof(int32_t), 48, 48, i32);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32_fast, sizeof(int32_t), 8, 8, i32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32_fast, sizeof(int32_t), 16, 16, i32);
    PROFILE_TRANSP(     1, isVerbose, mtx_transpose32x32_fast, sizeof(int32_t), 32, 32, i32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transpose32x32_fast, sizeof(int32_t), 48, 48, i32);
  }
}

static void mips_matop_phase2(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.i32;

  if (phaseNum == 0 || phaseNum == 2)
  {
    PROFILE_mpyf(     1, isVerbose, 16, 16, 16);
    PROFILE_mpyf(     1, isVerbose, 32, 32, 32);
    PROFILE_mpyf(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpyf(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpyf(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpyf(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpyf(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpyf(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpyf(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpyf(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpyf(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpyf_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpyf_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpyf_fast(isFull, isVerbose, 8, 80, 4);
    PROFILE_mpyf_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpyf_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpyf_fast(1, isVerbose, 8, 16, 4);
    PROFILE_mpyf_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mpytf(     1, isVerbose, 16, 16, 16);
    PROFILE_mpytf(     1, isVerbose, 32, 32, 32);
    PROFILE_mpytf(isFull, isVerbose, 40, 80, 8);
    PROFILE_mpytf(isFull, isVerbose, 40, 81, 8);
    PROFILE_mpytf(isFull, isVerbose, 40, 82, 8);
    PROFILE_mpytf(isFull, isVerbose, 40, 83, 8);
    PROFILE_mpytf(isFull, isVerbose, 2, 100, 8);
    PROFILE_mpytf(isFull, isVerbose, 8, 80, 2);
    PROFILE_mpytf(isFull, isVerbose, 8, 4, 2);
    PROFILE_mpytf(isFull, isVerbose, 8, 16, 2);
    PROFILE_mpytf(isFull, isVerbose, 8, 32, 2);
    PROFILE_mpytf_fast(1, isVerbose, 16, 16, 16);
    PROFILE_mpytf_fast(1, isVerbose, 32, 32, 32);
    PROFILE_mpytf_fast(isFull, isVerbose, 8, 80, 4);
    PROFILE_mpytf_fast(isFull, isVerbose, 8, 84, 4);
    PROFILE_mpytf_fast(isFull, isVerbose, 8, 4, 4);
    PROFILE_mpytf_fast(1, isVerbose, 8, 16, 4);
    PROFILE_mpytf_fast(isFull, isVerbose, 8, 32, 4);

    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf, 16, 100);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf, 16, 101);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf, 16, 102);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf, 16, 103);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf, 16, 104);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf, 40, 40);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf_fast, 16, 100);
    PROFILE_mtx_vecmpyf(1, isVerbose, mtx_vecmpyf_fast, 16, 104);
    PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf_fast, 40, 40);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef, sizeof(float32_t), 16, 16, f32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef, sizeof(float32_t), 27, 27, f32);
    PROFILE_TRANSP(     1, isVerbose, mtx_transposef, sizeof(float32_t), 32, 32, f32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef, sizeof(float32_t), 39, 39, f32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef, sizeof(float32_t), 48, 48, f32);

    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef_fast, sizeof(float32_t), 8, 8, f32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef_fast, sizeof(float32_t), 16, 16, f32);
    PROFILE_TRANSP(     1, isVerbose, mtx_transposef_fast, sizeof(float32_t), 32, 32, f32);
    PROFILE_TRANSP(isFull, isVerbose, mtx_transposef_fast, sizeof(float32_t), 48, 48, f32);
  }
}
void mips_matop(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {
      mips_matop_mpy(phaseNum, isFull, isVerbose, fout);
      mips_matop_vecmpy(phaseNum, isFull, isVerbose, fout);
      mips_matop_transp(phaseNum, isFull, isVerbose, fout);
    }
    if ( phaseNum == 0 || phaseNum == 2 )
    {
      mips_matop_phase2( phaseNum,  isFull,  isVerbose, fout);
    }
} /* mips_matop() */
