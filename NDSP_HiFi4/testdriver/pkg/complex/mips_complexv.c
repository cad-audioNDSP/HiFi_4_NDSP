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
* Test module for testing cycle performance (Vectorized Complex Math)
*/

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environemnt convfiguration. */
#include "config.h"
/* Library API */
#include LIBRARY_HEADER(complex)
/* Measurement utilities */
#include "mips.h"

void mips_complexv(int phaseNum, int isFull, int isVerbose, FILE * fout)
{

    if ( phaseNum == 0 || phaseNum == 1 )
    {
    }
    if ( phaseNum == 0 || phaseNum == 2 )
    {
        PROFILE_NORMALIZED(1, isVerbose, vec_complex2mag, (mips.out0.f32, mips.inp0.cf32, 200), fout, "N=200", prf_cyclespts, 200);
        PROFILE_NORMALIZED(1, isVerbose, vec_complex2invmag, (mips.out0.f32, mips.inp0.cf32, 200), fout, "N=200", prf_cyclespts, 200);
    }
} /* mips_complexv() */
