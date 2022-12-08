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
* Test module for testing cycle performance (Correlation, Convolution, 
* Dispreading, LMS)
*/

#include "config.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"


void mips_convol(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
     void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 60),fout,"N: 80; M: 60" ,prf_maccycle,    80*60 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 84),fout,"N: 256; M: 84",prf_maccycle,   84*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4,  80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16  ,  80, 60),fout,"N: 80; M: 60" ,prf_maccycle,    80*60 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4, 256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16  , 256, 84),fout,"N: 256; M: 84",prf_maccycle,   84*256 );
#if 0// for HiFi3/3z
    PROFILE_INVERTED(isFull, isVerbose,fir_convol24x24,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol24x24,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
#endif
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x32,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol32x32,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x32ep,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,  80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol32x32ep,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_convola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 56),fout,"N=80; M=56"   ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 80),fout,"N=256; M=80"  ,prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4,  80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4,  80, 60),fout,"N: 80; M: 60" ,prf_maccycle,    80*60 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4, 256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4, 256, 84),fout,"N: 256; M: 84",prf_maccycle,   84*256 );
#if 0// for HiFi3/3z
    PROFILE_INVERTED(isFull, isVerbose,fir_convola24x24,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola24x24,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
#endif
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N=80; M=56"   ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N=256; M=80"  ,prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x32ep,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,  80, 56),fout,"N=80; M=56"   ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola32x32ep,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80),fout,"N=256; M=80"  ,prf_maccycle,   80*256 );

}
void mips_cx_convol(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.u8;
  PROFILE_INVERTED(isFull, isVerbose,cxfir_convol32x16, (mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 2, 80, 56), fout, "N: 80; M: 56",prf_maccycle, 4*80 * 56);
  PROFILE_INVERTED(     1, isVerbose,cxfir_convol32x16, (mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 2, 256, 80), fout, "N: 256; M: 80",prf_maccycle, 4*80 * 256);
  PROFILE_INVERTED(isFull, isVerbose,cxfir_convola32x16, (pScr, mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 2, 80, 56), fout, "N: 80; M: 56",prf_maccycle, 4*80 * 56);
  PROFILE_INVERTED(     1, isVerbose,cxfir_convola32x16, (pScr, mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 2, 256, 80), fout, "N: 256; M: 80",prf_maccycle, 4*80 * 256);
}
void mips_xcorr(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.u8;
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorr16x16, (mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose, fir_xcorr16x16, (mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 60  ),fout, "N: 80; M: 60" , prf_maccycle, 80 * 60);
  PROFILE_INVERTED(     1, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 84 ),fout, "N: 256; M: 84", prf_maccycle, 84 * 256);
#if 0// for HiFi3/3z
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorr24x24, (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose, fir_xcorr24x24, (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
#endif
  PROFILE_INVERTED(isFull, isVerbose,fir_xcorr32x32,  (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose,fir_xcorr32x32,  (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
  PROFILE_INVERTED(isFull, isVerbose,fir_xcorr32x32ep,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose,fir_xcorr32x32ep,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
  PROFILE_INVERTED(isFull, isVerbose,cxfir_xcorr32x32,  (mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci32, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 4*80 * 56);
  PROFILE_INVERTED(     1, isVerbose,cxfir_xcorr32x32,  (mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci32, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 4*80 * 256);

  PROFILE_INVERTED(isFull, isVerbose,fir_xcorra16x16,  (pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle,    80*56);
  PROFILE_INVERTED(     1, isVerbose,fir_xcorra16x16,  (pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle,   80*256);
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 60  ),fout,"N: 80; M: 60" ,prf_maccycle, 80 * 60);
  PROFILE_INVERTED(     1, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 84 ),fout,"N: 256; M: 84",prf_maccycle, 84 * 256);
#if 0// for HiFi3/3z
  PROFILE_INVERTED(isFull, isVerbose, fir_xcorra24x24, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose, fir_xcorra24x24, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);
#endif
  PROFILE_INVERTED(isFull, isVerbose,fir_xcorra32x32,  (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose,fir_xcorra32x32,  (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);
  PROFILE_INVERTED(isFull, isVerbose,fir_xcorra32x32ep,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
  PROFILE_INVERTED(     1, isVerbose,fir_xcorra32x32ep,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);

}

void mips_blms(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp2.i16, mips.inp0.i16, 0x111, 111, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp2.i16, mips.inp0.i16, 0x111, 111, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );

    PROFILE_INVERTED(isFull, isVerbose, fir_blms16x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 80, 16), fout, "N: 80; M: 16", prf_maccycle, 80 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose, fir_blms16x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 64, 16), fout, "N: 64; M: 16", prf_maccycle, 64 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose, fir_blms16x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 64, 64), fout, "N: 64; M: 64", prf_maccycle, 64 * 64 * 2);
    PROFILE_INVERTED(isFull, isVerbose, fir_blms16x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 80, 64), fout, "N: 80; M: 64", prf_maccycle, 80 * 64 * 2);
    PROFILE_INVERTED(     1, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp2.i16, mips.inp0.i16, 0x111111, 111, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp2.i16, mips.inp0.i16, 0x111111, 111, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );
#if 0// for HiFi3/3z
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );
#endif
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 16), fout, "N: 80; M: 16", prf_maccycle, 80 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 16), fout, "N: 64; M: 16", prf_maccycle, 64 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 64), fout, "N: 64; M: 64", prf_maccycle, 64 * 64 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 64), fout, "N: 80; M: 64", prf_maccycle, 80 * 64 * 2);
    PROFILE_INVERTED(     1, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 80, 128), fout, "N: 80; M: 128", prf_maccycle, 80 * 128 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 64, 128), fout, "N: 64; M: 128", prf_maccycle, 64 * 128 * 2);

    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32ep, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 16), fout, "N: 80; M: 16", prf_maccycle, 80 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32ep, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 16), fout, "N: 64; M: 16", prf_maccycle, 64 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32ep, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 64), fout, "N: 64; M: 64", prf_maccycle, 64 * 64 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32ep, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 64), fout, "N: 80; M: 64", prf_maccycle, 80 * 64 * 2);
    PROFILE_INVERTED(     1, isVerbose,fir_blms32x32ep, (mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 80, 128), fout, "N: 80; M: 128", prf_maccycle, 80 * 128 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32ep, (mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 64, 128), fout, "N: 64; M: 128", prf_maccycle, 64 * 128 * 2);

    PROFILE_INVERTED(isFull, isVerbose,cxfir_blms32x32, (mips.out1.ci32, mips.out2.ci32, mips.inp1.ci32, mips.inp0.ci32, 0x111111, 111, 80, 16), fout, "N: 80; M: 16", prf_maccycle, 80 * 16 * 2*4);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blms32x32, (mips.out1.ci32, mips.out2.ci32, mips.inp1.ci32, mips.inp0.ci32, 0x111111, 111, 64, 16), fout, "N: 64; M: 16", prf_maccycle, 64 * 16 * 2*4);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blms32x32, (mips.out1.ci32, mips.out2.ci32, mips.inp1.ci32, mips.inp0.ci32, 0x111111, 111, 64, 64), fout, "N: 64; M: 64", prf_maccycle, 64 * 64 * 2*4);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blms32x32, (mips.out1.ci32, mips.out2.ci32, mips.inp1.ci32, mips.inp0.ci32, 0x111111, 111, 80, 64), fout, "N: 80; M: 64", prf_maccycle, 80 * 64 * 2*4);
    PROFILE_INVERTED(     1, isVerbose,cxfir_blms32x32, (mips.out1.ci32, mips.out2.ci32, mips.inp2.ci32, mips.inp0.ci32, 0x111111, 111, 80, 128), fout, "N: 80; M: 128", prf_maccycle, 80 * 128 * 2*4);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blms32x32, (mips.out1.ci32, mips.out2.ci32, mips.inp2.ci32, mips.inp0.ci32, 0x111111, 111, 64, 128), fout, "N: 64; M: 128", prf_maccycle, 64 * 128 * 2*4);

}

void mips_convolf(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_convolf,(mips.out0.f32, mips.inp2.f32, mips.inp1.f32,    80, 56),fout,"N: 80; M: 56",prf_maccycle,     80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convolf,(mips.out0.f32, mips.inp2.f32, mips.inp1.f32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_convolaf,(pScr, mips.out0.f32, mips.inp2.f32, mips.inp1.f32,    80, 56),fout,"N: 80; M: 56",prf_maccycle,     80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convolaf,(pScr, mips.out0.f32, mips.inp2.f32, mips.inp1.f32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
}

void mips_acorr(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr16x16,(mips.out0.i16, mips.inp2.i16, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr16x16,(mips.out0.i16, mips.inp2.i16, 256),fout,"N: 256",prf_maccycle,  256*256 );
#if 0// for HiFi3/3z
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr24x24,(mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr24x24,(mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );
#endif
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr32x32,(mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr32x32,(mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr32x32ep,(mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr32x32ep,(mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_acorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 80),fout,"N=80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 256),fout,"N=256",prf_maccycle,  256*256 );
#if 0// for HiFi3/3z
    PROFILE_INVERTED(isFull, isVerbose,fir_acorra24x24,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra24x24,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );
#endif
    PROFILE_INVERTED(isFull, isVerbose,fir_acorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorra32x32ep,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra32x32ep,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );

}

void mips_acorrf(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_acorrf,(mips.out0.f32, mips.inp2.f32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorrf,(mips.out0.f32, mips.inp2.f32, 256),fout,"N: 256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorraf,(pScr, mips.out0.f32, mips.inp2.f32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorraf,(pScr, mips.out0.f32, mips.inp2.f32, 256),fout,"N: 256",prf_maccycle,  256*256 );
}

void mips_xcorrf(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_xcorrf,  (mips.out0.f32, mips.inp2.f32, mips.inp1.f32, 80, 56    ),fout,"N: 80; M: 56" ,prf_maccycle,    80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorrf,  (mips.out0.f32, mips.inp2.f32, mips.inp1.f32, 256, 80   ),fout,"N: 256; M: 80",prf_maccycle,   80*256);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_xcorrf,(mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 80, 56 ),fout,"N: 80; M: 56" ,prf_maccycle,  4*80*56);
    PROFILE_INVERTED(     1, isVerbose,cxfir_xcorrf,(mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 256, 80),fout,"N: 256; M: 80",prf_maccycle, 4*80*256);

    PROFILE_INVERTED(isFull, isVerbose,fir_xcorraf,  (pScr, mips.out0.f32,  mips.inp2.f32,  mips.inp1.f32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle,    80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorraf,  (pScr, mips.out0.f32,  mips.inp2.f32,  mips.inp1.f32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle,   80*256);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_xcorraf,(pScr, mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 80, 56 ),fout,"N: 80; M: 56" ,prf_maccycle,  4*80*56);
    PROFILE_INVERTED(     1, isVerbose,cxfir_xcorraf,(pScr, mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 256, 80),fout,"N: 256; M: 80",prf_maccycle, 4*80*256);
}

void mips_blmsf(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp2.f32, mips.inp0.f32, 0.1f, 1.1f, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp2.f32, mips.inp0.f32, 0.1f, 1.1f, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );
}
void mips_cxblmsf(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2*4 );
    PROFILE_INVERTED(     1, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp2.cf32, mips.inp0.cf32, 0.1f, 1.1f, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp2.cf32, mips.inp0.cf32, 0.1f, 1.1f, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2*4 );
}

static void  mips_lacorr(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_lacorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 80),fout,"N=80",prf_maccycle,    80*80  /2 );
    PROFILE_INVERTED(     1, isVerbose,fir_lacorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 256),fout,"N=256",prf_maccycle,  256*256/2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_lacorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N=80",prf_maccycle,    80*80   /2);
    PROFILE_INVERTED(     1, isVerbose,fir_lacorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N=256",prf_maccycle,  256*256 /2);

}

static void  mips_lxcorr(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_lxcorra16x16,  (pScr,mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 80, 56    ),fout,"N=80; M=56" ,prf_maccycle,   80  * 56);
    PROFILE_INVERTED(     1, isVerbose,fir_lxcorra16x16,  (pScr,mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 256, 80   ),fout,"N=256; M=80",prf_maccycle,   256 * 80);
    PROFILE_INVERTED(isFull, isVerbose,fir_lxcorra32x32,  (pScr, mips.out0.i32,  mips.inp2.i32,  mips.inp1.i32, 80, 56  ),fout,"N=80; M=56" ,prf_maccycle,    80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_lxcorra32x32,  (pScr, mips.out0.i32,  mips.inp2.i32,  mips.inp1.i32, 256, 80 ),fout,"N=256; M=80",prf_maccycle,   80*256);

}

static void  mips_lconvol(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_lconvola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 56),fout,"N=80; M=56",prf_maccycle,     80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_lconvola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 80),fout,"N=256; M=80",prf_maccycle,   80*256);
    PROFILE_INVERTED(isFull, isVerbose,fir_lconvola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N=80; M=56",prf_maccycle,     80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_lconvola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N=256; M=80",prf_maccycle,   80*256);

}


void mips_firother(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
    if ( phaseNum == 0 || phaseNum == 1 )
    {
      mips_convol(phaseNum, isFull, isVerbose, fout);
      mips_cx_convol(phaseNum, isFull, isVerbose, fout);
      mips_lconvol(phaseNum, isFull, isVerbose, fout);
      mips_xcorr(phaseNum, isFull, isVerbose, fout);
      mips_lxcorr(phaseNum, isFull, isVerbose, fout);
      mips_acorr(phaseNum, isFull, isVerbose, fout);
      mips_lacorr(phaseNum, isFull, isVerbose, fout);
      mips_blms(phaseNum, isFull, isVerbose, fout);
    }
    if ( phaseNum == 0 || phaseNum == 2 )
    {
      mips_convolf(phaseNum, isFull, isVerbose, fout);
      mips_xcorrf(phaseNum, isFull, isVerbose, fout);
      mips_acorrf(phaseNum, isFull, isVerbose, fout);
      mips_blmsf(phaseNum, isFull, isVerbose, fout);
      mips_cxblmsf(phaseNum, isFull, isVerbose, fout);
    }
} /* mips_firother() */
