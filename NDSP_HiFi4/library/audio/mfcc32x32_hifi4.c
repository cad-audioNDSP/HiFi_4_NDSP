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
  NatureDSP Signal Processing Library. Audio processing part
    Compute Mel-Frequency Cepstrum Coefficients (MFCC) from speech signal
    32-bit fixed-point variant
    C code optimized for HiFi4 core
  IntegrIT, 2006-2018
*/

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_audio.h"
#include "NatureDSP_Signal_matop.h"
/* MFCC extractor internal declarations. */
#include "mfcc_internal.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)
#define ALIGN_PAD      (ALIGN_SIZE-1)
#define ALIGN_PTR(p)   (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)
#define sz_i32         sizeof(fract32)

/*-------------------------------------------------------------------------
  Compute Mel-Frequency Cepstrum Coefficients (MFCC) from speech signal.
  MFCC features extraction procedure comprises the following steps:
  - input speech signal is optionally passed through a pre-emphasis filter;
  - filtered signal is subject to short-time Fourier transform (STFT) followed by
    magnitude spectrum computation;
  - a set of filters is applied to the magnitude spectrum, with triangular weight
    functions constructed in such a way that the prescribed frequency range 
    is divided into overlapping bands of equal mel-frequency width;
  - log-scaled filterbank energies are decorrelated via a Discrete Cosine Transform 
    Type II (DCT-II) to form cepstrum coefficients;
  - in the last step a sine lifter is optionally applied to cepstra to align 
    coefficient magnitudes.
  In general, the computation procedure follows the MFCC features extraction 
  algorithm adopted in the Hidden Markov Models Toolkit (HTK), as descibed in:
  [1] S. Young, G. Evermann, M. Gales, T. Hain, D. Kershaw, X. Liu, G. Moore,
      J. Odell, D. Ollason, D. Povey, V. Valtchev, P. Woodland,
      The HTK Book (for HTK version 3.4), 
      Cambridge University Engineering Department, 2009. 
      http://htk.eng.cam.ac.uk/docs/docs.shtml
  In addition, a number of options provide an ability to emulate the operation of
  another popular package for speech analysis:
  [2] The Auditory Toolbox for MATLAB by Malcolm Slaney, Version 2
      Interval Research Corporation
      https://engineering.purdue.edu/~malcolm/interval/1998-010/
  Precision: 
  32x32                       32-bit fixed-point input/output data
  f                           Single precision floating-point input/output data
  Input:
  objmem                      Memory block allocated for the instance object:
  params                      MFCC features extraction parameters
  callback                    User-supplied callback functions
  speech[stftHopLen]          Speech samples; Q31 for 32x32
  Temporary:
  pScr                        Scratch memory area for the processing function. To 
                              determine the scratch area size, use the respective
                              helper function: mfcc<32x32|f>_getScratchSize()
  Output:
  cepstra[cepstraNum]         Cepstral coefficients, the number of fractional bits
                              for 32x32 is defined by MFCC_CEPSTRA_FRACT_BITS.
  Restrictions:
  cepstra[],speech[]          Must not overlap, and must be aligned by 8-bytes
  Fs                          8000 <= Fs <= 48000
  fftSize                     256, 512, 1024 or 2048
  stftWinLen, stftHopLen      Must be multiples of 2, 0 < stftHopLen <= stftWinLen <= fftSize
  mfbLowFreqQ8, mfbUppFreqQ8  0 <= mfbLowFreqQ8 < mfbUppFreqQ8 <= 16000*256
  mfbBandNum, cepstraNum      0 < cepstraNum <= mfbBandNum <= 40
-------------------------------------------------------------------------*/

/* Returns: size of memory to be allocated for an instance object, in bytes. */
size_t mfcc32x32_alloc( const mfcc_params_t * params )
{
    logmel_params_t logmel_params;
    int dctRows, dctCols;
    size_t allocSize, allocSize_logmel;
    NASSERT(NULL!=params);
    /* Check if mfcc parameters meet the restrictions stated by the API. */
    if (!(8000<=params->Fs && params->Fs<=48000)                                                                   ||
        !(params->fftSize==256 || params->fftSize==512 || params->fftSize==1024 || params->fftSize==2048)          ||
        !((params->stftHopLen%2)==0 && 0<params->stftHopLen && params->stftHopLen<=params->stftWinLen)             ||
        !((params->stftWinLen%2)==0 && params->stftWinLen<=params->fftSize)                                        ||
        !(0<=params->mfbLowFreqQ8 && params->mfbLowFreqQ8<params->mfbUppFreqQ8 && params->mfbUppFreqQ8<=16000*256) ||
        !(0<params->cepstraNum && params->cepstraNum<=params->mfbBandNum)                                          ||
        !(params->mfbBandNum<=40) 
       ) {
        NASSERT(!"mfcc32x32_alloc(): bad params"); return (0);
    }
    /* Fast matrix-vector multiply routine requires that matrix dimensions are divisible by 4. */
    dctRows = MFCC32X32_ALWAYS_USE_MTXVECMPY_FAST ? (params->cepstraNum+3)&~3 : params->cepstraNum;
    dctCols = MFCC32X32_ALWAYS_USE_MTXVECMPY_FAST ? (params->mfbBandNum+3)&~3 : params->mfbBandNum;
    /* Compute the allocation size. */
    allocSize = ALIGN_PAD + sizeof(mfcc32x32_t)       + /* Instance structure */
                ALIGN_PAD + sz_i32*params->stftWinLen + /* stftWeights        */
                ALIGN_PAD + sz_i32*params->stftWinLen + /* stftFrame          */
                ALIGN_PAD + sz_i32*dctRows*dctCols;     /* dctMtx             */
    if (params->lifter>0) {
        allocSize += ALIGN_PAD + sz_i32*params->cepstraNum;  /* lifterCoefs   */
    }
    LOGMEL_PARAMS_FROM_MFCC_PARAMS(&logmel_params, params);
    LOGMEL_PARAMS_FROM_MFCC_PARAMS(&logmel_params, params);
    if (0!=(allocSize_logmel = logmel32x32_alloc(&logmel_params))) {
        return allocSize+allocSize_logmel;
    } else {
        return 0;
    }
} /* mfcc32x32_alloc() */

/* Returns: handle to the object, or NULL if initialization failed. */
mfcc32x32_handle_t mfcc32x32_init( void * objmem, const mfcc_params_t * params, const mfcc32x32_callback_t * callback )
{
    void * p = objmem, *logmel_objmem;
    mfcc32x32_t * mfcc;
    fract32 *stftWeights, *stftFrame, *dctMtx, *lifterCoefs=NULL;
    size_t allocSize;
    logmel_params_t logmel_params;
    logmel32x32_handle_t logmel;
    size_t logmel_allocSize;
    int dctRows, dctCols, dctExp, lifterExp=0;
    int n;
    NASSERT(NULL!=objmem && NULL!=params && NULL!=callback);
    /* Validate MFCC parameters and recalculate the instance object size. */
    allocSize = mfcc32x32_alloc(params); (void)allocSize;
    LOGMEL_PARAMS_FROM_MFCC_PARAMS(&logmel_params, params);
    logmel_allocSize = logmel32x32_alloc(&logmel_params);
    /* DCT matrix dimensions. */
    dctRows = MFCC32X32_ALWAYS_USE_MTXVECMPY_FAST ? (params->cepstraNum+3)&~3 : params->cepstraNum;
    dctCols = MFCC32X32_ALWAYS_USE_MTXVECMPY_FAST ? (params->mfbBandNum+3)&~3 : params->mfbBandNum;
    /* Partition the object memory space. */
    mfcc        = (mfcc32x32_t*)ALIGN_PTR(p); p = mfcc+1;
    stftWeights = (fract32    *)ALIGN_PTR(p); p = stftWeights + params->stftWinLen;
    stftFrame   = (fract32    *)ALIGN_PTR(p); p = stftFrame + params->stftWinLen;
    dctMtx      = (fract32    *)ALIGN_PTR(p); p = dctMtx + dctRows*dctCols;
    if (params->lifter>0) {
        lifterCoefs = (fract32*)ALIGN_PTR(p); p = lifterCoefs + params->cepstraNum;
    }
    logmel_objmem = p; p = (uint8_t*)logmel_objmem + logmel_allocSize;
    NASSERT((int8_t*)p-(int8_t*)objmem <= (int)allocSize);
    /* Marginal initialization for the memory debugger. */
    AE_S32X2_X(AE_ZERO32(), (ae_int32x2*)stftFrame, params->stftWinLen*sz_i32);
    /* Make an internal copy of window weights, */
    if (NULL!=callback->genWindow) {
        callback->genWindow(callback->host, stftWeights, params->stftWinLen);
    } else {
        stftWeights = NULL;
    }
    /* Wipe the STFT signal frame. */
    for ( n=0; n<params->stftWinLen; n++ ) stftFrame[n] = 0;
    logmel = logmel32x32_init(logmel_objmem, &logmel_params);
    if (NULL==logmel) return NULL;
    /* Compute the DCT matrix. */
    dctExp = mfcc32x32_compDctMatrix(dctMtx, params->cepstraNum, params->mfbBandNum, 
                                     dctRows, dctCols,  MFCC_OPT_EQ(params->opt, DCT, ORTHOGONAL));
    /* Optionally compute the lifter coefficients. */
    if (params->lifter>0) {
        lifterExp = mfcc32x32_compLifterCoefs(lifterCoefs, params->lifter, params->cepstraNum);
    }
    /* Initialize the instance object. */
    memset(mfcc, 0, sizeof(*mfcc));
    mfcc->magic        = MFCC32X32_MAGIC;
    mfcc->params       = *params;
    mfcc->callback     = *callback;
    mfcc->logmel       = logmel;
    mfcc->stftWeights  = stftWeights;
    mfcc->stftFrame    = stftFrame;
    mfcc->preemphState = 0;
    mfcc->dctRows      = dctRows;
    mfcc->dctCols      = dctCols;
    mfcc->dctMtx       = dctMtx;
    mfcc->dctExp       = dctExp;
    mfcc->lifterCoefs  = lifterCoefs;
    mfcc->lifterExp    = lifterExp;
    /* Select either the fast or regular variant of the matrix-vector multiply routine. */
    mfcc->mtxvecmpy = (0==((dctRows|dctCols)&3)) ? &mtx_vecmpy32x32_fast : &mtx_vecmpy32x32;
    return (mfcc32x32_handle_t)mfcc;
} /* mfcc32x32_init() */
