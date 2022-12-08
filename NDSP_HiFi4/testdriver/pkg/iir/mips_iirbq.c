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
* Test module for testing cycle performance (Biquad Filters)
*/

#include "config.h"
#include LIBRARY_HEADER(iir)
#include "mips.h"


/* IIR performance measurement tests */
//----bqriir---
static const int16_t coef_sos_16[5*8] =
{ // b0      b1      b2      a1     a2
   16384,  31383,  16384,  12682, 14622,
   16384, -31383,  16384, -12682, 14622,
   16384,  32215,  16384,   7625, 11691,
   16384, -32215,  16384,  -7625, 11691,
   16384,      0, -16384,      0, 10377,
   16384,  31383,  16384,  12682, 14622,
   16384, -31383,  16384, -12682, 14622,
   16384,  32215,  16384,   7625, 11691,
};
static const int32_t coef_sos_32[5*8] =
{ //  b0          b1         b2         a1         a2
  1073741824, 2056704919, 1073741824, 831104644,958261518,
  1073741824,-2056704919, 1073741824,-831104644,958261518,
  1073741824, 2111239901, 1073741824, 499713750,766176384,
  1073741824,-2111239901, 1073741824,-499713750,766176384,
  1073741824,          0,-1073741824,         0,680063938,
  1073741824, 2056704919, 1073741824, 831104644,958261518,
  1073741824,-2056704919, 1073741824,-831104644,958261518,
  1073741824, 2111239901, 1073741824, 499713750,766176384,
};
static const float32_t coef_sos_f[5*16] =
{ //  b0          b1         b2         a1         a2
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f,
1.0000000f,-1.9662454f, 1.0000000f,-0.4653947f, 0.7135574f,
1.0000000f, 0.0000000f,-1.0000000f, 0.0000000f, 0.6333589f,
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f,
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f,
1.0000000f,-1.9662454f, 1.0000000f,-0.4653947f, 0.7135574f,
1.0000000f, 0.0000000f,-1.0000000f, 0.0000000f, 0.6333589f,
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f
};

static const int16_t coef_g[8] =
{
  2460,19682,9107,9107,22170,2460,19682,9107
};


#define OBJ_PROFILE_NORMALIZED_IIR(_cond,_verb,_objname, _a_arg, _i_arg, _p_arg, _file, _info_,_fmt, _norm)  \
{                                                                        \
   _objname##_handle_t handle=NULL;                                      \
  int isPresent;                                                         \
  isPresent =IS_PRESENT(_objname##_alloc);               \
  isPresent|=IS_PRESENT(_objname##_init);                \
  isPresent |= IS_PRESENT(_objname);                     \
  if (isPresent )     handle = _objname##_init _i_arg;                   \
  if (handle == NULL) isPresent = 0;                                     \
  PROFILE_NORMALIZED(_cond,_verb,_objname, _p_arg, _file, _info_, _fmt, _norm)       \
}

/* IIR performance measurement tests */
#define PROFILE_BQRIIR16X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir16x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR16X16_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir16x16_df2_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X16_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x16_df2_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR24x24_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir24x24_df1_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR24x24_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir24x24_df2_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X32_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x32_df1_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X32_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x32_df2_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF1_ND(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df1_nd ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF2_ND(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df2_nd ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF2T_ND(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df2t_nd,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQCIIRF_DF1_ND(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqciirf_df1_nd ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.cf32, mips.inp1.cf32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR16X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir16x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR32X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir32x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR32X32_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir32x32_df1_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIRF_DF1_ND(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriirf_df1_nd ,(M),(objinstance_memory,M,coef_sos_f,1,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);

#define PROFILE_BQRIIR16X16_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir16x16_df1,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR16X16_DF2(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir16x16_df2,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X16_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x16_df1,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X16_DF2(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x16_df2,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR24x24_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir24x24_df1,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR24x24_DF2(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir24x24_df2,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X32_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x32_df1,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X32_DF2(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x32_df2,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF1(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df1 ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF2(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df2 ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF2T(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df2t,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQCIIRF_DF1(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqciirf_df1 ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.cf32, mips.inp1.cf32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR16X16_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir16x16_df1,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR32X16_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir32x16_df1,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR32X32_DF1(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir32x32_df1,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIRF_DF1(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriirf_df1 ,(M),(objinstance_memory,M,coef_sos_f,1,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);

#define PROFILE_NEWIIR(fun,N,M,iirstate_type,suffixc,suffixd)                                           \
{                                                                                                       \
    iirstate_type state;                                                                                \
    iir_init((&state),mips.inp2.suffixc,mips.inp1.suffixd,M);                                                     \
    PROFILE_NORMALIZED(fun,(&state,mips.out0.suffixd,mips.inp0.suffixd,N),fout,"N=" #N ",M=" #M,prf_cyclesbqd,(N*M)); \
}

#define PROFILE_IIRDF1f(fun,N,M)     PROFILE_NEWIIR(fun,N,M,iir_statef,f32)

void mips_iirbq(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR16X16_DF1(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR16X16_DF1(isFull, isVerbose, 80, 5, 1)
                                     
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR16X16_DF2(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR16X16_DF2(isFull, isVerbose, 80, 5, 1)

        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 1, 0)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 2, 1)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 3, 0)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 4, 1)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 5, 0)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 6, 1)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,256, 7, 0)
        PROFILE_BQRIIR32X16_DF1(     1, isVerbose,256, 8, 1)
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,80, 5, 0) 
        PROFILE_BQRIIR32X16_DF1(isFull, isVerbose,80, 5, 1)
   
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR32X16_DF2(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR32X16_DF2(isFull, isVerbose, 80, 5, 1)
#if 0 //HiFi3/3z API
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR24x24_DF1(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR24x24_DF1(isFull, isVerbose, 80, 5, 1)
                                     
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR24x24_DF2(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR24x24_DF2(isFull, isVerbose, 80, 5, 1)
#endif
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR32X32_DF1(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR32X32_DF1(isFull, isVerbose, 80, 5, 1)
    
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 1, 0) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 2, 1) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 3, 0) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 4, 1) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 5, 0) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 6, 1) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 256, 7, 0) 
        PROFILE_BQRIIR32X32_DF2(     1, isVerbose, 256, 8, 1) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 80, 5, 0) 
        PROFILE_BQRIIR32X32_DF2(isFull, isVerbose, 80, 5, 1)

        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 1, 0) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 2, 1) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 3, 0) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 4, 1) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 5, 0) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 6, 1) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 256, 7, 0) 
        PROFILE_STEREO_BQRIIR16X16_DF1(     1, isVerbose, 256, 8, 1) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 80, 5, 0) 
        PROFILE_STEREO_BQRIIR16X16_DF1(isFull, isVerbose, 80, 5, 1)

        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 1, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 2, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 3, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 4, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 5, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 6, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,256, 7, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1(     1, isVerbose,256, 8, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,80, 5, 0) 
        PROFILE_STEREO_BQRIIR32X16_DF1(isFull, isVerbose,80, 5, 1)

        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 1, 0) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 2, 1) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 3, 0) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 4, 1) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 5, 0) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 6, 1) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 256, 7, 0) 
        PROFILE_STEREO_BQRIIR32X32_DF1(     1, isVerbose, 256, 8, 1) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 80, 5, 0) 
        PROFILE_STEREO_BQRIIR32X32_DF1(isFull, isVerbose, 80, 5, 1)

    }

    if ( phaseNum == 0 || phaseNum == 2 )
    {
        PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 1) 
        PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 2) 
        PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 3) 
        PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 4) 
        PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 8) 
        PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512,12) 
        PROFILE_BQRIIRF_DF1(     1, isVerbose, 512,16) 

        PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 1) 
        PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 2) 
        PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 3) 
        PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 4) 
        PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 8) 
        PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512,12) 
        PROFILE_BQRIIRF_DF2(     1, isVerbose, 512,16) 

        PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 1) 
        PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 2) 
        PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 3) 
        PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 4) 
        PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 8) 
        PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512,12) 
        PROFILE_BQRIIRF_DF2T(     1, isVerbose, 512,16) 

        PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 1) 
        PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 2) 
        PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 3) 
        PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 4) 
        PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 8) 
        PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512,12) 
        PROFILE_BQCIIRF_DF1(     1, isVerbose, 512,16) 

        PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 1) 
        PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 2) 
        PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 3) 
        PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 4) 
        PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 8) 
        PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512,12) 
        PROFILE_STEREO_BQRIIRF_DF1(     1, isVerbose, 512,16) 

    }
} /* mips_iirbq() */


void mips_iirbq_nd(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_BQRIIR16X16_DF1_ND(1, isVerbose, 256, 8, 1)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_BQRIIR16X16_DF2_ND(1, isVerbose, 256, 8, 1)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_BQRIIR32X16_DF1_ND(1, isVerbose, 256, 8, 1)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_BQRIIR32X16_DF2_ND(1, isVerbose, 256, 8, 1)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_BQRIIR32X32_DF1_ND(1, isVerbose, 256, 8, 1)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_BQRIIR32X32_DF2_ND(1, isVerbose, 256, 8, 1)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(1, isVerbose, 256, 8, 1)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(1, isVerbose, 256, 8, 1)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 1, 0)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 2, 1)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 3, 0)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 4, 1)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 5, 0)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 6, 1)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 7, 0)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(1, isVerbose, 256, 8, 1)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 0)
        PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 1)
    }

    if ( phaseNum == 0 || phaseNum == 2 )
    {
        PROFILE_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 1) 
        PROFILE_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 2)
        PROFILE_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 3)
        PROFILE_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 4)
        PROFILE_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 8)
        PROFILE_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 12)
        PROFILE_BQRIIRF_DF1_ND(1, isVerbose, 512, 16)

        PROFILE_BQRIIRF_DF2_ND(isFull, isVerbose, 512, 1)
        PROFILE_BQRIIRF_DF2_ND(isFull, isVerbose, 512, 2)
        PROFILE_BQRIIRF_DF2_ND(isFull, isVerbose, 512, 3)
        PROFILE_BQRIIRF_DF2_ND(isFull, isVerbose, 512, 4)
        PROFILE_BQRIIRF_DF2_ND(isFull, isVerbose, 512, 8)
        PROFILE_BQRIIRF_DF2_ND(isFull, isVerbose, 512, 12)
        PROFILE_BQRIIRF_DF2_ND(1, isVerbose, 512, 16)

        PROFILE_BQRIIRF_DF2T_ND(isFull, isVerbose, 512, 1)
        PROFILE_BQRIIRF_DF2T_ND(isFull, isVerbose, 512, 2)
        PROFILE_BQRIIRF_DF2T_ND(isFull, isVerbose, 512, 3)
        PROFILE_BQRIIRF_DF2T_ND(isFull, isVerbose, 512, 4)
        PROFILE_BQRIIRF_DF2T_ND(isFull, isVerbose, 512, 8)
        PROFILE_BQRIIRF_DF2T_ND(isFull, isVerbose, 512, 12)
        PROFILE_BQRIIRF_DF2T_ND(1, isVerbose, 512, 16)

        PROFILE_BQCIIRF_DF1_ND(isFull, isVerbose, 512, 1)
        PROFILE_BQCIIRF_DF1_ND(isFull, isVerbose, 512, 2)
        PROFILE_BQCIIRF_DF1_ND(isFull, isVerbose, 512, 3)
        PROFILE_BQCIIRF_DF1_ND(isFull, isVerbose, 512, 4)
        PROFILE_BQCIIRF_DF1_ND(isFull, isVerbose, 512, 8)
        PROFILE_BQCIIRF_DF1_ND(isFull, isVerbose, 512, 12)
        PROFILE_BQCIIRF_DF1_ND(1, isVerbose, 512, 16)

        PROFILE_STEREO_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 1)
        PROFILE_STEREO_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 2)
        PROFILE_STEREO_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 3)
        PROFILE_STEREO_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 4)
        PROFILE_STEREO_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 8)
        PROFILE_STEREO_BQRIIRF_DF1_ND(isFull, isVerbose, 512, 12)
        PROFILE_STEREO_BQRIIRF_DF1_ND(1, isVerbose, 512, 16)
    }
} /* mips_iirbq_nd() */
