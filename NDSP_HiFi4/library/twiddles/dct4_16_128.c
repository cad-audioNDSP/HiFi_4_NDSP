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
    NatureDSP Signal Processing Library. FFT part
    DCT-IV twiddles
    Integrit, 2006-2018
*/

#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "dct4_twd.h"

static const complex_fract16 ALIGN(8) dct4_twd_128[]=
{
    {{32767,-201  }},{{32762,-603  }},{{32753,-1005 }},{{32738,-1407 }},
    {{32718,-1809 }},{{32693,-2210 }},{{32664,-2611 }},{{32629,-3012 }},
    {{32590,-3412 }},{{32546,-3812 }},{{32496,-4211 }},{{32442,-4609 }},
    {{32383,-5007 }},{{32319,-5404 }},{{32251,-5800 }},{{32177,-6195 }},
    {{32099,-6590 }},{{32015,-6983 }},{{31927,-7376 }},{{31834,-7767 }},
    {{31737,-8157 }},{{31634,-8546 }},{{31527,-8933 }},{{31415,-9319 }},
    {{31298,-9704 }},{{31177,-10088}},{{31050,-10469}},{{30920,-10850}},
    {{30784,-11228}},{{30644,-11605}},{{30499,-11980}},{{30350,-12354}},
    {{23028,-23312}},{{22740,-23593}},{{22449,-23870}},{{22154,-24144}},
    {{21856,-24414}},{{21555,-24680}},{{21251,-24943}},{{20943,-25202}},
    {{20632,-25457}},{{20318,-25708}},{{20001,-25956}},{{19681,-26199}},
    {{19358,-26439}},{{19032,-26674}},{{18703,-26906}},{{18372,-27133}},
    {{18037,-27357}},{{17700,-27576}},{{17361,-27791}},{{17018,-28002}},
    {{16673,-28209}},{{16326,-28411}},{{15976,-28610}},{{15624,-28803}},
    {{15269,-28993}},{{14912,-29178}},{{14553,-29359}},{{14192,-29535}},
    {{13828,-29707}},{{13463,-29875}},{{13095,-30038}},{{12725,-30196}}
};

static const complex_fract16 dct3_128[]=
{
    {{32767,0    }},{{32758,804  }},{{32729,1608 }},{{32679,2411 }},
    {{32610,3212 }},{{32522,4011 }},{{32413,4808 }},{{32286,5602 }},
    {{32138,6393 }},{{31972,7180 }},{{31786,7962 }},{{31581,8740 }},
    {{31357,9512 }},{{31114,10279}},{{30853,11039}},{{30572,11793}},
    {{30274,12540}},{{29957,13279}},{{29622,14010}},{{29269,14733}},
    {{28899,15447}},{{28511,16151}},{{28106,16846}},{{27684,17531}},
    {{27246,18205}},{{26791,18868}},{{26320,19520}},{{25833,20160}},
    {{25330,20788}},{{24812,21403}},{{24279,22006}},{{23732,22595}}
};

static const complex_fract16 rfft_64[]=
{
    {{32610,3212 }},
    {{32138,6393 }},
    {{31357,9512 }},
    {{30274,12540}},
    {{28899,15447}},
    {{27246,18205}},
    {{25330,20788}},
    {{23170,23170}},
    {{20788,25330}},
    {{18205,27246}},
    {{15447,28899}},
    {{12540,30274}},
    {{9512 ,31357}},
    {{6393 ,32138}},
    {{3212 ,32610}}
};

static const complex_fract16 fft_32[]=
{
    {{ 32767, 0    }},
    {{ 32767, 0    }},
    {{ 32767, 0    }},
    {{ 32138,-6393 }},
    {{ 30274,-12540}},
    {{ 27246,-18205}},
    {{ 30274,-12540}},
    {{ 23170,-23170}},
    {{ 12540,-30274}},
    {{ 27246,-18205}},
    {{ 12540,-30274}},
    {{-6393 ,-32138}},
    {{ 23170,-23170}},
    {{ 0    ,-32768}},
    {{-23170,-23170}},
    {{ 18205,-27246}},
    {{-12540,-30274}},
    {{-32138,-6393 }},
    {{ 12540,-30274}},
    {{-23170,-23170}},
    {{-30274, 12540}},
    {{ 6393 ,-32138}},
    {{-30274,-12540}},
    {{-18205, 27246}}
};

static const tdct4_twd_fr16 descr = {128,dct4_twd_128,dct3_128,rfft_64 ,fft_32 };
const dct_handle_t dct4_16_128=(dct_handle_t)&descr;
const dct_handle_t mdct_16_128=(dct_handle_t)&descr;
