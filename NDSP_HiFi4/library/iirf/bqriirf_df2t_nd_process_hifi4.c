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
  NatureDSP Signal Processing Library. IIR part
    Bi-quad Real Block IIR, floating point, Direct Form II transposed
    C code optimized for HiFi4
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Bi-quad Real Block IIR
  Computes a real IIR filter (cascaded IIR direct form I or II using 5 
  coefficients per bi-quad + gain term). Real input data are stored
  in vector x. The filter output result is stored in vector r. The filter 
  calculates N output samples using SOS and G matrices.
  NOTE:
  1. Bi-quad coefficients may be derived from standard SOS and G matrices 
  generated by MATLAB. However, typically biquad stages have big peaks 
  in their step response which may cause undesirable overflows at the 
  intermediate outputs. To avoid that the additional scale factors 
  coef_g[M] may be applied. These per-section scale factors may require 
  some tuning to find a compromise between quantization noise and possible 
  overflows. Output of the last section is directed to an additional 
  multiplier, with the gain factor being a power of two, either negative 
  or non-negative. It is specified through the total gain shift amount 
  parameter gain of each filter initialization function.
  2. 16x16 filters may suffer more from accumulation of the roundoff errors,
  so filters should be properly designed to match noise requirements
  3. Due to the performance reasons, IIR biquad filters may introduce 
  additional algorithmic delay of several sampless. Amount of that delay
  might be requested by the  xxx_groupDelay API. For sensitive applications
  all the filters have delayless implementations (with  _nd  suffix in the name).
  Formally, the xxx_groupDelay APIs is also implemented for that kind of filters,
  but return zero.
  
  Precision: 
  16x16         16-bit data, 16-bit coefficients, 16-bit intermediate 
                stage outputs (DF1, DF1 stereo, DF II form)
  32x16         32-bit data, 16-bit coefficients, 32-bit intermediate 
                (DF1, DF1 stereo, DF II form) stage outputs
  32x32         32-bit data, 32-bit coefficients, 32-bit intermediate 
                (DF I, DF1 stereo,  DF II form) stage outputs 
  f             floating point (DF I, DF1 stereo, DF II and DF IIt)

  Input:
  N             length of input sample block
  M             number of bi-quad sections
  S             1 for mono, 2 for stereo API
  s[]           scratch memory area (for fixed-point functions only). 
                Minimum number of bytes depends on selected filter structure 
                and precision:
                  BQRIIR16X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR16X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF2_SCRATCH_SIZE( N,  M ),
                  STEREO_BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIRF_DF1_SCRATCH_SIZE    ( N, M )
                 If a particular macro returns zero, then the corresponding
                 IIR doesn't require a scratch area and parameter s may 
                 hold zero.
  coef_sos[M*5]  filter coefficients stored in blocks of 5 numbers: 
                 b0 b1 b2 a1 a2. 
                 For fixed-point funcions, fixed point format of filter 
                 coefficients is Q1.14 for 16x16 and 32x16, or Q1.30 for 32x32 
  coef_sosl[M*5] filter coefficients for the left channel (stereo filters only)
  coef_sosr[M*5] filter coefficients for the left channel (stereo filters only)
  coef_g[M]      scale factor for each section, Q15 (for fixed-point 
                 functions only)
  coef_gl[M]     scale factor for the left channel (stereo filters only)
  coef_gr[M]     scale factor for the right channel (stereo filters only)
  gain           total gain shift amount applied to output signal of the
                 last section, -48..15
  gainl          total gain shift amount  for the left channel (stereo filters 
                 only)
  gainr          total gain shift amount for the right channel (stereo filters 
                 only)

  x[N*S]         input samples, Q15, Q31 or floating point. Stereo samples 
                 go in interleaved order (left, right)
  Output:
  r[N*S]         output data, Q15, Q31 or floating point. Stereo samples go 
                 in interleaved order (left, right) 

  Restriction:
  x,r,s,coef_g,coef_sos should not overlap
  N   - must be a multiple of 2
  s[] - whenever supplied must be aligned on an 8-bytes boundary
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bqriirf_df2t_nd.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,bqriirf_df2t_nd,( bqriirf_df2t_nd_handle_t _bqriir,
                      float32_t * restrict       z,
                const float32_t *            x, int N))
#elif (HAVE_VFPU)

// Process data. The filter instance is identified with a handle returned by
// the initializing function.
void bqriirf_df2t_nd(bqriirf_df2t_nd_handle_t _bqriir,
    float32_t *  restrict       z,
    const float32_t *                 x, int N)
{
    bqriirf_df2t_nd_t *state;
    const ae_int32x2 * restrict pX;
    ae_int32x2 * restrict pZ;
    const xtfloatx2  * restrict pb;
    const xtfloat  * restrict pa1,
        *restrict pa2,
        *restrict pb0,
        *restrict pb1,
        *restrict pb2;
    const xtfloatx2  * restrict pDrd;
    xtfloatx2  * restrict pDwr;
    //xtfloatx2 t0;
    xtfloatx2 dx0, dz0, dw0, t0;
    xtfloatx2 dx1, dz1, dw1, t1;
    xtfloatx2 a10, a20, b00, b10, b20;    
    xtfloatx2 a11, a21, b01, b11, b21;
    xtfloatx2 tz;
    ae_int32x2 tmp;
    int n, m;
    int M;
    NASSERT(_bqriir);
    if (N <= 0) return;
    NASSERT(N % 2 == 0);
    NASSERT(x);
    NASSERT(z);
    state = (bqriirf_df2t_nd_t*)(_bqriir);
    NASSERT(state);
    NASSERT(state->st);
    NASSERT(state->cf);
    M = state->M;

    /* Initialize pointers */
    pX = (const ae_int32x2 *)(x);
    pb = (const xtfloatx2  *)(state->cf);
    pDrd = (const xtfloatx2  *)(state->st);
    pDwr = (xtfloatx2  *)(state->st);

    /* Process sections */
    for (m = 0; m < (M>>2); m++)
    {
        pZ = (ae_int32x2 *)z;
        /* load delay lines and coefficients */
        XT_LSX2IP(dz0, pDrd, 2*sizeof(float32_t));
        XT_LSX2IP(dw0, pDrd, 2*sizeof(float32_t));
        XT_LSX2IP(dz1, pDrd, 2 * sizeof(float32_t));
        XT_LSX2IP(dw1, pDrd, 2 * sizeof(float32_t));
        /* load part of the coefficients */
        b00 = XT_LSX2I(pb, 0 * sizeof(float32_t));
        b01 = XT_LSX2I(pb, 2 * sizeof(float32_t));
        b10 = XT_LSX2I(pb, 4 * sizeof(float32_t));
        b11 = XT_LSX2I(pb, 6 * sizeof(float32_t));
        
        /* processing loop */
        for (n = 0; n<(N); n++)
        {           
            /* load filter coefficients */
            b20 = XT_LSX2I(pb, 8 * sizeof(float32_t));
            b21 = XT_LSX2I(pb, 10 * sizeof(float32_t));
            a10 = XT_LSX2I(pb, 12 * sizeof(float32_t));
            a11 = XT_LSX2I(pb, 14 * sizeof(float32_t));
            a20 = XT_LSX2X(pb, 16 * sizeof(float32_t));
            a21 = XT_LSX2X(pb, 18 * sizeof(float32_t));

            tz = AE_ZERO32();
            /* load input sample */
            AE_L32_IP(tmp, castxcc(ae_int32, pX), sizeof(float32_t));
            dx0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
            t0 = dz0;
            XT_MADDN_SX2(t0, dx0, b00);
            dx0 = XT_SEL32_HH_SX2(dx0, t0);
            t0 = dz0;
            XT_MADDN_SX2(t0, dx0, b00);
            dz0 = dw0;
            XT_MADDN_SX2(dz0, b10, dx0);
            XT_MSUBN_SX2(dz0, a10, t0);
            dw0 = XT_MUL_SX2(b20, dx0);
            XT_MSUBN_SX2(dw0, a20, t0);

            dx1 = XT_SEL32_LH_SX2(t0, tz);
            t1 = dz1;
            XT_MADDN_SX2(t1, dx1, b01);
            dx1 = XT_SEL32_HH_SX2(dx1, t1);
            t1 = dz1;
            XT_MADDN_SX2(t1, dx1, b01);
            dz1 = dw1;
            XT_MADDN_SX2(dz1, b11, dx1);
            XT_MSUBN_SX2(dz1, a11, t1);
            dw1 = XT_MUL_SX2(b21, dx1);
            XT_MSUBN_SX2(dw1, a21, t1);

            /* save output */
            tmp = XT_AE_MOVINT32X2_FROMXTFLOATX2(t1);
            AE_S32_L_IP(tmp, castxcc(ae_int32, pZ), sizeof(float32_t));          
        }
        pb += 10;
        XT_SSX2IP(dz0, pDwr, 2*sizeof(float32_t));
        XT_SSX2IP(dw0, pDwr, 2*sizeof(float32_t));
        XT_SSX2IP(dz1, pDwr, 2 * sizeof(float32_t));
        XT_SSX2IP(dw1, pDwr, 2 * sizeof(float32_t));
        /* switch pointer to the input data to the pointer to the output */
        pX = (const ae_int32x2 *)(z);
    }
    if (M & 2)
    {
        for (m = 0; m < 1; m++)
        {
            pZ = (ae_int32x2 *)z;
            /* load delay lines and coefficients */
            XT_LSX2IP(dz0, pDrd, 2 * sizeof(float32_t));
            XT_LSX2IP(dw0, pDrd, 2 * sizeof(float32_t));
            /* load filter coefficients */
            b00 = XT_LSX2I(pb, 0 * sizeof(float32_t));           
            b10 = XT_LSX2I(pb, 2 * sizeof(float32_t));  
            b20 = XT_LSX2I(pb, 4 * sizeof(float32_t));
            a10 = XT_LSX2I(pb, 6 * sizeof(float32_t));
            a20 = XT_LSX2I(pb, 8 * sizeof(float32_t));
            /* processing loop */
            for (n = 0; n<(N); n++)
            {                              
                tz = AE_ZERO32();
                /* load input sample */
                AE_L32_IP(tmp, castxcc(ae_int32, pX), sizeof(float32_t));
                dx0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
                //dx0 = XT_SEL32_HH_SX2(dx0, tz);
                t0 = dz0;
                XT_MADDN_SX2(t0, dx0, b00);
                dx0 = XT_SEL32_HH_SX2(dx0, t0);
                t0 = dz0;
                XT_MADDN_SX2(t0, dx0, b00);
                dz0 = dw0;
                XT_MADDN_SX2(dz0, b10, dx0);
                XT_MSUBN_SX2(dz0, a10, t0);
                dw0 = XT_MUL_SX2(b20, dx0);
                XT_MSUBN_SX2(dw0, a20, t0);
                /* save output */
                tmp = XT_AE_MOVINT32X2_FROMXTFLOATX2(t0);
                AE_S32_L_IP(tmp, castxcc(ae_int32, pZ), sizeof(float32_t));
            }
            pb += 5;
            XT_SSX2IP(dz0, pDwr, 2 * sizeof(float32_t));
            XT_SSX2IP(dw0, pDwr, 2 * sizeof(float32_t));
            /* switch pointer to the input data to the pointer to the output */
            pX = (const ae_int32x2 *)(z);
        }
    }
    if (M & 1)
    {
        pb0 = (const xtfloat *)((float32_t *)pb);
        pb1 = (const xtfloat *)((float32_t *)pb0 + 1);
        pb2 = (const xtfloat *)((float32_t *)pb1 + 1);
        pa1 = (const xtfloat *)((float32_t *)pb2 + 1);
        pa2 = (const xtfloat *)((float32_t *)pa1 + 1);

        /* Process sections */
        for (m = 0; m < M%2; m++)
        {
            xtfloat dx, dz, dw, t;
            xtfloat a1, a2, b0, b1, b2;
            pZ = (ae_int32x2 *)z;
            /* load delay lines and coefficients */
            XT_LSIP(dz, castxcc(xtfloat, pDrd), sizeof(float32_t));
            XT_LSIP(dw, castxcc(xtfloat, pDrd), sizeof(float32_t));
            XT_LSIP(b0, pb0, 5 * sizeof(float32_t));
            XT_LSIP(b1, pb1, 5 * sizeof(float32_t));
            XT_LSIP(b2, pb2, 5 * sizeof(float32_t));
            XT_LSIP(a1, pa1, 5 * sizeof(float32_t));
            XT_LSIP(a2, pa2, 5 * sizeof(float32_t));

            /* processing loop */
            for (n = 0; n<(N >> 1); n++)
            {
                XT_LSIP(dx, castxcc(xtfloat, pX), sizeof(float32_t));
                t = dz;
                XT_MADDN_S(t, dx, b0);
                dz = dw;
                XT_MADDN_S(dz, b1, dx);
                XT_MSUBN_S(dz, a1, t);
                dw = XT_MUL_S(b2, dx);
                XT_MSUBN_S(dw, a2, t);
                XT_SSIP(t, castxcc(xtfloat, pZ), sizeof(float32_t));
                XT_LSIP(dx, castxcc(xtfloat, pX), sizeof(float32_t));
                t = dz;
                XT_MADDN_S(t, b0, dx);
                dz = dw;
                XT_MADDN_S(dz, b1, dx);
                XT_MSUBN_S(dz, a1, t);
                dw = XT_MUL_S(b2, dx);
                XT_MSUBN_S(dw, a2, t);
                XT_SSIP(t, castxcc(xtfloat, pZ), sizeof(float32_t));
            }

            XT_SSIP(dz, castxcc(xtfloat, pDwr), sizeof(float32_t));
            XT_SSIP(dw, castxcc(xtfloat, pDwr), sizeof(float32_t));
            /* switch pointer to the input data to the pointer to the output */
            pX = (const ae_int32x2 *)(z);
        }
    }


    /* final scaling */
    {
        xtfloatx2 ft;
        int32_t s;
        s = state->gain;
        s = ((s + 127) & 255) << 23;
        ft = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVDA32X2(s, s));
        pZ = (ae_int32x2 *)z;

        if (((uintptr_t)z)&(sizeof(ae_int32x2)-1))
        {
            xtfloat t0, ft0;
            ft0 = XT_WFR(s);
            /* output pointer is not aligned by 8-byte boundary */
            XT_LSIP(t0, castxcc(xtfloat, pX), sizeof(xtfloat));
            t0 = XT_MUL_S(t0, ft0);
            XT_SSIP(t0, castxcc(xtfloat, pZ), sizeof(xtfloat));
            t0 = XT_LSX((xtfloat *)pX, (N - 2)*sizeof(xtfloat));
            t0 = XT_MUL_S(t0, ft0);
            XT_SSX(t0, (xtfloat *)pZ, (N - 2)*sizeof(xtfloat));
            N -= 2;
        }
        for (n = 0; n < (N >> 1); n++)
        {
            XT_LSX2IP(t0, castxcc(xtfloatx2, pX), sizeof(xtfloatx2));
            t0 = t0*ft;
            XT_SSX2IP(t0, castxcc(xtfloatx2, pZ), sizeof(xtfloatx2));
        }
    }
}/* bqriirf_df2t_nd_process() */
#else 
// code for scalar FPU
void bqriirf_df2t_nd( bqriirf_df2t_nd_handle_t _bqriir,
                   float32_t *  restrict       z,
             const float32_t *                 x, int N)
{
    bqriirf_df2t_nd_t *state;
    const xtfloat *cf;
    const xtfloat * restrict pDrd;
          xtfloat * restrict pDwr;
    const xtfloat * pX;
          xtfloat * pZ;
    const xtfloat * restrict pXr;
          xtfloat * restrict pZr;
    float32_t scale;
    int32_t s;
    int n,m;
    int M;
    NASSERT(_bqriir);
    state=(bqriirf_df2t_nd_t*)(_bqriir);
    if(N<=0) return;
    NASSERT(N%2==0);
    NASSERT(x);
    NASSERT(z);
    NASSERT(state);
    NASSERT(state->st);
    NASSERT(state->cf);
    M=state->M;
    cf=(const xtfloat*)state->cf;
    pDrd=(const xtfloat*)state->st;
    pDwr=(      xtfloat*)state->st;
    s=state->gain;
    s=((s+127)&255)<<23;
    scale=XT_WFR(s);

    for (m=0; m<M; m++,x=z)
    {
        xtfloat a1,a2,b0,b1,b2,dx,dz,dw;
        pX=(const xtfloat*)x;
        pZ=(      xtfloat*)z;
        XT_LSIP(dz,pDrd,sizeof(xtfloat));
        XT_LSIP(dw,pDrd,sizeof(xtfloat));
        XT_LSIP(b0,cf,sizeof(xtfloat));
        XT_LSIP(b1,cf,sizeof(xtfloat));
        XT_LSIP(b2,cf,sizeof(xtfloat));
        XT_LSIP(a1,cf,sizeof(xtfloat));
        XT_LSIP(a2,cf,sizeof(xtfloat));
        for (n=0;n<N;n+=2)
        {
            xtfloat t,u0,u1;
            XT_LSIP(u0,pX,sizeof(xtfloat));
            XT_LSIP(u1,pX,sizeof(xtfloat));
            
            t=dz;  XT_MADD_S(t,b0,u0);
            dz=dw; XT_MADD_S(dz,b1,u0);
            XT_MSUB_S(dz,a1,t) ;
            dw =XT_MUL_S(b2,u0);
            XT_MSUB_S(dw,a2,t);
            dx =t;
            XT_SSIP(dx, pZ, sizeof(xtfloat));
            
            t=dz;  XT_MADD_S(t,b0,u1);
            dz=dw; XT_MADD_S(dz,b1,u1);
            XT_MSUB_S(dz,a1,t) ;
            dw =XT_MUL_S(b2,u1);
            XT_MSUB_S(dw,a2,t);
            dx =t;
            XT_SSIP(dx, pZ, sizeof(xtfloat));
        }
        XT_SSIP(dz,pDwr,sizeof(xtfloat));
        XT_SSIP(dw,pDwr,sizeof(xtfloat));
    }
    __Pragma("no_reorder")
    /* final scaling */
    {
        pZr=(      xtfloat*)z;
        pXr=(const xtfloat*)x;
        for (n=0; n<(N>>1); n++)
        {
            xtfloat x0,x1;
            XT_LSIP(x0,pXr,sizeof(xtfloat));
            XT_LSIP(x1,pXr,sizeof(xtfloat));
            XT_SSIP(XT_MUL_S(x0,scale),pZr,sizeof(xtfloat));
            XT_SSIP(XT_MUL_S(x1,scale),pZr,sizeof(xtfloat));
        }
    }
}

#endif