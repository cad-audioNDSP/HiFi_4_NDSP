/* ------------------------------------------------------------------------ */
/* Copyright (c) 2019 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2015-2019 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
* Test module for testing cycle performance (image processing APIs.)
*/

#include <string.h>
#include <math.h>
/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(img)
#include "NatureDSP_Signal_img.h"
/* MIPS measurement means. */
#include "mips.h"
/* Utility functions and macros. */
#include "utils.h"
#include <stdlib.h>
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define TEST_ODD_ANGLES 1

static const int32_t rgbyuv_transform[]={160524403,315143225,61203284,-78989817,-155080532,234075718,329638740,-276483151,-53692460,611941572,-211876105,-311755570,1090980749};

#define PROFILE_IMG_FFT(isFull,isVerbose,basename,_file,resolution,w,h) { \
        size_t szScr;                                                     \
        size_t bytesOut;                                                  \
        void* pScr;                                                       \
        void* pImg;                                                       \
        void* pOutImg;                                                       \
        int istride = (w+7)&~7;                                           \
        imgsize_t sz;                                                     \
        sz.width=w;                                                       \
        sz.height=h;                                                      \
        sz.stride=istride;                                                \
        szScr=basename##_getScratchSize (&sz );                           \
        szScr=(szScr+15)&~15;                                             \
        /* bytesOut=h*w*sizeof(complex_fract16); */                            \
        bytesOut=h*w*sizeof(fract16);                             \
        NASSERT(szScr+3 * bytesOut<=sizeof(mips));                            \
        (void)(bytesOut);                                                  \
        pScr=(void*)&mips;                                                \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                            \
        pOutImg=(void*)(((uintptr_t)pImg)+bytesOut);                            \
        PROFILE_SIMPLE(isFull,isVerbose,basename,                         \
                      /* (pScr,(complex_fract16*)pImg,pImg,128<<16,&sz), */   \
                       (pScr,(complex_fract16*)pOutImg,pImg,128<<16,&sz),    \
                       _file,resolution,prf_cycle);                       \
}

#define PROFILE_IMG_IFFT(isFull,isVerbose,basename,_file,resolution,w,h) {\
        size_t szScr;                                                     \
        size_t bytesIn;                                                   \
        void* pScr;                                                       \
        void* pImg;                                                       \
        void* pOutImg;                                                       \
        int istride = (w+7)&~7;                                           \
        imgsize_t sz;                                                     \
        sz.width=w;                                                       \
        sz.height=h;                                                      \
        sz.stride=istride;                                                \
        szScr=basename##_getScratchSize (&sz );                           \
        szScr=(szScr+15)&~15;                                             \
        /*bytesIn=h*w*sizeof(complex_fract16);  */                            \
        bytesIn=h*w*sizeof(fract16);                              \
        (void)(bytesIn);                                                  \
        NASSERT(szScr+3 * bytesIn<=sizeof(mips));                             \
        pScr=(void*)&mips;                                                \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                            \
        pOutImg=(void*)(((uintptr_t)pImg)+2 * bytesIn);                            \
        PROFILE_SIMPLE(isFull,isVerbose,basename,                         \
                      /* (pScr, pImg, (complex_fract16*)pImg,128<<16,&sz,0), */ \
                       (pScr, pOutImg, (complex_fract16*)pImg,128<<16,&sz,0),\
                       _file,resolution,prf_cycle);                       \
}


#define PROFILE_IMG_RESIZE(isFull,isVerbose,basename,fast,method,_file,resolution,win,hin,wout,hout,method_str,szpixel) {  \
        size_t szObj, szScr;                                         \
        void* memObj;                                                \
        void* memScr;                                                \
        int istride = fast?(win+7)&~7:win+1;                         \
        int ostride = fast?(wout+7)&~7:wout+1;                       \
        imgresize_params_t params={{win,hin,istride},                \
                                   {wout,hout,ostride},method};      \
        imgresize_handle_t handle;                                   \
        szObj=basename##_alloc  (&params);                           \
        szScr=basename##_getScratchSize (&params );                  \
        NASSERT(hin*istride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2)); \
        NASSERT(hout*ostride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));\
        /*NASSERT(szScr<sizeof(tProfiler_scratch)+sizeof(mips.inp0)); */ \
        (void)szScr;                                                   \
        memObj=malloc(szObj);                                          \
        memScr=mallocAlign(szScr, 4 * NatureDSP_Signal_get_isa_opt(NATUREDSP_ISA_OPT_INT16_SIMD_WIDTH));  \
        handle=basename##_init(memObj,&params);                        \
        PROFILE_SIMPLE(isFull,isVerbose,basename##_process,            \
                       (handle,memScr,(void*)&mips.out0,\
                       (const void*)&mips.inp0),                       \
                       _file,resolution ", " method_str,prf_cycle);    \
        freeAlign(memScr);                                                  \
        free(memObj);                                                  \
}
  
#define PROFILE_IMG_ROTATE(isFull,isVerbose,basename,fast,_file,resolution,w,h,angle,szpixel) {  \
        size_t szObj, szScr;                                                \
        imgsize_t szOut;                                                    \
        int bytesIn,bytesOut;                                               \
        void* memObj;                                                       \
        void* pScr;                                                         \
        void* pImg;                                                         \
        int stride = fast?((w+7)&~7):w+1;                                   \
        imgrotate_params_t params;                                          \
        imgrotate_handle_t handle;                                          \
        params.in.width=w;                                                  \
        params.in.height=h;                                                 \
        params.in.stride=stride;                                            \
        params.fill=0;                                                      \
        params.angleQ15=(int16_t)(((int64_t)11930464L*angle+32768)>>16);    \
        szObj=basename##_alloc  (&params);                                  \
        szScr=basename##_getScratchSize (&params );                         \
        basename##_getOutSize (&szOut,&params );                            \
        szScr=(szScr+15)&~15;                                               \
        bytesIn=szpixel*w*stride;                                           \
        bytesOut=szOut.width*szOut.stride;                                  \
        NASSERT(szScr + MAX(bytesIn,bytesOut)<sizeof(mips));                  \
        (void)szObj,(void)szScr,(void)bytesIn,(void)bytesOut;               \
        pScr=(void*)&mips;                                                  \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                              \
        memObj=(void*)objinstance_memory;                                   \
        NASSERT(szObj<sizeof(objinstance_memory));                          \
        handle=basename##_init(memObj,&params);                             \
        PROFILE_SIMPLE(isFull,isVerbose,basename##_process,                 \
                       (handle,pScr,pImg,pImg,&szOut),                      \
                       _file,resolution " " #angle " degrees",prf_cycle);   \
}

#define PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull,isVerbose,basename,fast,_file,resolution,w,h,angle,szpixel) {  \
        size_t szObj, szScr;                                                \
        imgsize_t szOut;                                                    \
        int bytesIn,bytesOut;                                               \
        void* memObj;                                                       \
        void* pScr;                                                         \
        void* pImg;                                                         \
        void* pOutImg;                                                         \
        int stride = fast?((w+7)&~7):w+1;                                   \
        imgrotate_params_t params;                                          \
        imgrotate_handle_t handle;                                          \
        params.in.width=w;                                                  \
        params.in.height=h;                                                 \
        params.in.stride=stride;                                            \
        params.fill=0;                                                      \
        params.angleQ15=(int16_t)(((int64_t)11930464L*angle+32768)>>16);    \
        szObj=basename##_alloc  (&params);                                  \
        szScr=basename##_getScratchSize (&params );                         \
        basename##_getOutSize (&szOut,&params );                            \
        szScr=(szScr+15)&~15;                                               \
        bytesIn=szpixel*w*stride;                                           \
        bytesOut=szOut.width*szOut.stride;                                  \
        NASSERT(szScr + (bytesIn + bytesOut)<sizeof(mips));                  \
        (void)szObj,(void)szScr,(void)bytesIn,(void)bytesOut;               \
        pScr=(void*)&mips;                                                  \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                              \
        pOutImg=(void*)(((uintptr_t)pImg)+bytesIn);                              \
        memObj=(void*)objinstance_memory;                                   \
        NASSERT(szObj<sizeof(objinstance_memory));                          \
        handle=basename##_init(memObj,&params);                             \
        PROFILE_SIMPLE(isFull,isVerbose,basename##_process,                 \
                       /*(handle,pScr,pImg,pImg,&szOut), */                     \
                        (handle,pScr,pOutImg,pImg,&szOut),                      \
                       _file,resolution " " #angle " degrees",prf_cycle);   \
}

#define PROFILE_IMG_HIST(isFull,isVerbose,f_name,fast,_file,resolution,w,h) {    \
        int stride = fast?(w+7)&~7:w+1;                                          \
        imgsize_t sz={w,h,stride};                                               \
        imghist_t hist={0,0,0,mips.out1.i32};                                    \
        NASSERT(h*stride<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2)); \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,(&hist,                           \
                       (void*)mips.inp0.u8,&sz,256),_file,resolution,prf_cycle); \
        }

#define PROFILE_IMG_NORM(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) {   \
        int stride = fast?(w+7)&~7:w+1;                                                 \
        imgsize_t sz={w,h,stride};                                                      \
        NASSERT(h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));\
        NASSERT(h*stride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));\
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                    \
                       (void*)mips.inp0.u8,&sz,0,255),_file,resolution,prf_cycle);      \
        }
#define PROFILE_IMG_NORM_NONLINEAR(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) {\
        int stride = fast?(w+7)&~7:w+1;                                                        \
        imgsize_t sz={w,h,stride};                                                             \
        NASSERT(h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));       \
        NASSERT(h*stride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));       \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                           \
                       (void*)mips.inp0.u8,&sz,mips.inp1.i16),_file,resolution,prf_cycle);     \
        }

#define PROFILE_IMG_INTERLEAVE(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) { \
        int stride = fast?(w+7)&~7:w+1;                                                     \
        imgsize_t sz={w,h,stride};                                                          \
        NASSERT(3*h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+   \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));  \
        NASSERT(  h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+   \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));  \
        /* PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.inp0.u8,    */                    \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                        \
                                                (const void*)mips.inp0.u8,                  \
                                                (const void*)mips.inp0.u8,                  \
                                                (const void*)mips.inp0.u8,&sz),             \
                                               _file,resolution,prf_cycle);                 \
        }
#define PROFILE_IMG_DEINTERLEAVE(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) {\
        int stride = fast?(w+7)&~7:w+1;                                                      \
        imgsize_t sz={w,h,stride};                                                           \
        NASSERT(3*h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+    \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));   \
        NASSERT(  h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+    \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));   \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                         \
                                                (void*)mips.out0.u8,                         \
                                                (void*)mips.out0.u8,                         \
                                                (const void*)mips.inp0.u8, &sz),             \
                                               _file,resolution,prf_cycle);                  \
        }

#define PROFILE_IMG_CONVERT(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) { \
        int stride = fast?(w+7)&~7:w+1;                                                  \
        imgsize_t sz={w,h,stride};                                                       \
        NASSERT(h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)); \
        NASSERT(h*stride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2)); \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                     \
                                                (void*)mips.out0.u8,                     \
                                                (void*)mips.out0.u8,                     \
                                                (const void*)mips.inp0.u8,               \
                                                (const void*)mips.inp1.u8,               \
                                                (const void*)mips.inp2.u8,               \
                                                rgbyuv_transform,&sz),                   \
                                               _file,resolution,prf_cycle);              \
        }

#define PROFILE_IMG_PAD(isFull,isVerbose,f_name,fast,_file,resolution,win,hin,wout,hout,szpixel) \
{                                                                           \
        int istride = fast?(win+7)&~7:win+1;                                \
        int ostride = fast?(wout+7)&~7:wout+1;                              \
        imgsize_t szin ={win,hin,istride};                                  \
        imgsize_t szout={wout,hout,ostride};                                \
        int x=(win-wout)>>1;                                                \
        int y=(hin-hout)>>1;                                                \
        imgpad_params_t params={szin,szout,x,y,-1};                         \
        size_t sz=f_name##_getScratchSize(&params);                         \
        NASSERT(sz<=sizeof(mips.scratch0));                                 \
        NASSERT(hout*ostride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));\
        NASSERT(hin*istride *szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));\
        (void)sz;                                                           \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.scratch0.u8,    \
                                                (void*)mips.out0.u8,        \
                                                (const void*)mips.inp0.u8,  \
                                                &params),                   \
                                               _file,resolution,prf_cycle); \
}


static void mips_hist(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"SQCIF(128x96)",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"QCIF(176x144)",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"CIF (352x288)",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"QVGA(320x240)",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"VGA (640x480)",640,480);

    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"SQCIF(128x96)",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"QCIF(176x144)",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"CIF (352x288)",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"QVGA(320x240)",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"VGA (640x480)",640,480);

    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "SQCIF(128x96)", 128, 96);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "QCIF(176x144)", 176, 144);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "CIF (352x288)", 352, 288);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "QVGA(320x240)", 320, 240);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "VGA (640x480)", 640, 480);

    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "SQCIF(128x96)", 128, 96);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "QCIF(176x144)", 176, 144);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "CIF (352x288)", 352, 288);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "QVGA(320x240)", 320, 240);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "VGA (640x480)", 640, 480);

    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"SQCIF(128x96)",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"QCIF(176x144)",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"CIF (352x288)",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"QVGA(320x240)",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"VGA (640x480)",640,480);

    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"SQCIF(128x96)",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"QCIF(176x144)",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"CIF (352x288)",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"QVGA(320x240)",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"VGA (640x480)",640,480);
}

static void mips_interleave(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imginterleave  ,0,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imginterleave16,0,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"QVGA(320x240)",320,240,2);
//    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imgfastinterleave  ,1,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imgfastinterleave16,1,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"QVGA(320x240)",320,240,2);
//    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgdeinterleave  ,0,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgdeinterleave16,0,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"QVGA(320x240)",320,240,2);
//    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgfastdeinterleave  ,1,fout,"CIF(352x288)",352,288,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgfastdeinterleave16,1,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"QVGA(320x240)",320,240,2);
//    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"VGA (640x480)",640,480,2);
}

static void mips_convert(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_rgbyuv  ,0,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_rgbyuv16,0,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_rgbyuv  ,1,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_rgbyuv16,1,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_yuvrgb  ,0,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_yuvrgb16,0,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_yuvrgb  ,1,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"VGA (640x480)",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_yuvrgb16,1,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"VGA (640x480)",640,480,2);
}

static void mips_norm(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"VGA (640x480)",640,480,1);

    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"VGA (640x480)",640,480,1);

    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "QCIF(176x144)", 176, 144, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "CIF (352x288)", 352, 288, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "QVGA(320x240)", 320, 240, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "VGA (640x480)", 640, 480, 1);

    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "QCIF(176x144)", 176, 144, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "CIF (352x288)", 352, 288, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "QVGA(320x240)", 320, 240, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "VGA (640x480)", 640, 480, 1); 

    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"VGA (640x480)",640,480,2);
    // nonlinear
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"VGA (640x480)",640,480,1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"SQCIF(128x96)",128,96,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"QCIF(176x144)",176,144,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"CIF (352x288)",352,288,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"QVGA(320x240)",320,240,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"VGA (640x480)",640,480,1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "SQCIF(128x96)", 128, 96, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "QCIF(176x144)", 176, 144, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "CIF (352x288)", 352, 288, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "QVGA(320x240)", 320, 240, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "VGA (640x480)", 640, 480, 1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "SQCIF(128x96)", 128, 96, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "QCIF(176x144)", 176, 144, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "CIF (352x288)", 352, 288, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "QVGA(320x240)", 320, 240, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "VGA (640x480)", 640, 480, 1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"VGA (640x480)",640,480,2);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"SQCIF(128x96)",128,96,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"QCIF(176x144)",176,144,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"CIF (352x288)",352,288,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"QVGA(320x240)",320,240,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"VGA (640x480)",640,480,2);
}

static void mips_imgrotate_gu8(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (0 == phaseNum || 1 == phaseNum)
  {
    (void)isFull;
    //8-bit unsigned
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "SQCIF(128x96)", 128, 96, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "QCIF(176x144)", 176, 144, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "CIF(352x288)", 352, 288, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gu8, 0, fout, "QVGA(320x240)", 320, 240, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "VGA(640x480)", 640, 480, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gu8, 0, fout, "VGA(640x480)", 640, 480, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "VGA(640x480)", 640, 480, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gu8, 0, fout, "VGA(640x480)", 640, 480, 270, 1);
#if TEST_ODD_ANGLES
    //    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gu8,0,fout,"VGA(640x480)",640,480, 45,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA(640x480)",640,480,135,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA(640x480)",640,480,225,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA(640x480)",640,480,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "SQCIF(128x96)", 128, 96, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QCIF(176x144)", 176, 144, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "CIF(352x288)", 352, 288, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "QVGA(320x240)", 320, 240, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "VGA(640x480)", 640, 480, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gu8, 1, fout, "VGA(640x480)", 640, 480, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "VGA(640x480)", 640, 480, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gu8, 1, fout, "VGA(640x480)", 640, 480, 270, 1);
#if TEST_ODD_ANGLES
    //    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gu8,1,fout,"VGA(640x480)",640,480, 45,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA(640x480)",640,480,135,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA(640x480)",640,480,225,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA(640x480)",640,480,315,1);
#endif

  }
} /* mips_imgrotate_gu8() */

static void mips_imgrotate_gs8(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (0 == phaseNum || 1 == phaseNum)
  {
    (void)isFull;
    //8-bit signed
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "SQCIF(128x96)", 128, 96, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "QCIF(176x144)", 176, 144, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "CIF(352x288)", 352, 288, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs8, 0, fout, "QVGA(320x240)", 320, 240, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "VGA(640x480)", 640, 480, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs8, 0, fout, "VGA(640x480)", 640, 480, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "VGA(640x480)", 640, 480, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs8, 0, fout, "VGA(640x480)", 640, 480, 270, 1);
#if TEST_ODD_ANGLES
    //    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs8,0,fout,"VGA(640x480)",640,480, 45,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA(640x480)",640,480,135,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA(640x480)",640,480,225,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA(640x480)",640,480,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "SQCIF(128x96)", 128, 96, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QCIF(176x144)", 176, 144, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "CIF(352x288)", 352, 288, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 270, 1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 45, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 135, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 225, 1);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "QVGA(320x240)", 320, 240, 315, 1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "VGA(640x480)", 640, 480, 0, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs8, 1, fout, "VGA(640x480)", 640, 480, 90, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "VGA(640x480)", 640, 480, 180, 1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs8, 1, fout, "VGA(640x480)", 640, 480, 270, 1);
#if TEST_ODD_ANGLES
    //    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs8,1,fout,"VGA(640x480)",640,480, 45,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA(640x480)",640,480,135,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA(640x480)",640,480,225,1);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA(640x480)",640,480,315,1);
#endif

  }
} /* mips_imgrotate_gs8() */

static void mips_imgrotate_gs16(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (0 == phaseNum || 1 == phaseNum)
  {
    (void)isFull;
    // 16-bit
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "SQCIF(128x96)", 128, 96, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "QCIF(176x144)", 176, 144, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "CIF(352x288)", 352, 288, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgrotate_gs16, 0, fout, "QVGA(320x240)", 320, 240, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "VGA(640x480)", 640, 480, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgrotate_gs16, 0, fout, "VGA(640x480)", 640, 480, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "VGA(640x480)", 640, 480, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgrotate_gs16, 0, fout, "VGA(640x480)", 640, 480, 270, 2);
#if TEST_ODD_ANGLES
    //    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs16,0,fout,"VGA(640x480)",640,480, 45,2);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA(640x480)",640,480,135,2);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA(640x480)",640,480,225,2);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA(640x480)",640,480,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "SQCIF(128x96)", 128, 96, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QCIF(176x144)", 176, 144, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "CIF(352x288)", 352, 288, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 270, 2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(     1, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 45, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 135, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 225, 2);
    PROFILE_IMG_ROTATE(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "QVGA(320x240)", 320, 240, 315, 2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "VGA(640x480)", 640, 480, 0, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(     1, isVerbose, imgfastrotate_gs16, 1, fout, "VGA(640x480)", 640, 480, 90, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "VGA(640x480)", 640, 480, 180, 2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose, imgfastrotate_gs16, 1, fout, "VGA(640x480)", 640, 480, 270, 2);
#if TEST_ODD_ANGLES
    //    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs16,1,fout,"VGA(640x480)",640,480, 45,2);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA(640x480)",640,480,135,2);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA(640x480)",640,480,225,2);
    //    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA(640x480)",640,480,315,2);
#endif

  }
} /* mips_imgrotate_gs16() */

void mips_imgrotate(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  mips_imgrotate_gu8(phaseNum, isFull, isVerbose, fout);
  mips_imgrotate_gs8(phaseNum, isFull, isVerbose, fout);
  mips_imgrotate_gs16(phaseNum, isFull, isVerbose, fout);
}

static void mips_resize_nearest_gu8(int isFull, int isVerbose, FILE * fout)
{
  // 8-bit unsigned
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_nearest, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "nearest", 1);

  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_nearest, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "nearest", 1);

}

static void mips_resize_nearest_gs8(int isFull, int isVerbose, FILE * fout)
{

  // 8-bit signed
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "nearest", 1);

  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "nearest", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "nearest", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "nearest", 1);

}

static void mips_resize_nearest_gs16(int isFull, int isVerbose, FILE * fout)
{
    // 16-bit 
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"nearest",2);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"nearest",2);

}

static void mips_resize_nearest(int isFull, int isVerbose, FILE * fout)
{
  mips_resize_nearest_gu8(isFull, isVerbose, fout);
  mips_resize_nearest_gs8(isFull, isVerbose, fout);
  mips_resize_nearest_gs16(isFull, isVerbose, fout);
}

static void mips_resize_bilinear_gu8(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gu8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "bilinear", 1);

  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gu8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "bilinear", 1);

}

static void mips_resize_bilinear_gs8(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_bilinear, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "bilinear", 1);

  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "bilinear", 1);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "bilinear", 1);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_bilinear, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "bilinear", 1);

}

static void mips_resize_bilinear_gs16(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs16, 0, imgresize_method_bilinear, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "bilinear", 2);

  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->QCIF(176x144)", 128, 96, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->CIF (352x288)", 128, 96, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->QVGA(320x240)", 128, 96, 320, 240, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "SQCIF(128x96)->VGA (640x480)", 128, 96, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->SQCIF(128x96)", 176, 144, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->CIF (352x288)", 176, 144, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->QVGA(320x240)", 176, 144, 320, 240, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QCIF(176x144)->VGA (640x480)", 176, 144, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "CIF (352x288)->SQCIF(128x96)", 352, 288, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "CIF (352x288)->QCIF(176x144)", 352, 288, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "CIF (352x288)->QVGA(320x240)", 352, 288, 320, 240, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "CIF (352x288)->VGA (640x480)", 352, 288, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->SQCIF(128x96)", 320, 240, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->QCIF(176x144)", 320, 240, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->CIF (352x288)", 320, 240, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "QVGA(320x240)->VGA (640x480)", 320, 240, 640, 480, "bilinear", 2);
  PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "VGA (640x480)->SQCIF(128x96)", 640, 480, 128, 96, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "VGA (640x480)->QCIF(176x144)", 640, 480, 176, 144, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "VGA (640x480)->CIF (352x288)", 640, 480, 352, 288, "bilinear", 2);
  PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs16, 1, imgresize_method_bilinear, fout, "VGA (640x480)->QVGA(320x240)", 640, 480, 320, 240, "bilinear", 2);

}

static void mips_resize_bilinear(int isFull, int isVerbose, FILE * fout)
{
  mips_resize_bilinear_gu8(isFull, isVerbose, fout);
  mips_resize_bilinear_gs8(isFull, isVerbose, fout);
  mips_resize_bilinear_gs16(isFull, isVerbose, fout);
}

static void mips_resize_bicubic_gu8(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"bicubic",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"bicubic",1);

}

static void mips_resize_bicubic_gs8(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,    isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"bicubic",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"bicubic",1);

}

static void mips_resize_bicubic_gs16(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"bicubic",2);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->QCIF(176x144)",128,96, 176,144, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->CIF (352x288)",128,96, 352,288, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->QVGA(320x240)",128,96, 320,240, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF(128x96)->VGA (640x480)",128,96, 640,480, "bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF(176x144)->SQCIF(128x96)",176,144, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF(176x144)->CIF (352x288)",176,144, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF(176x144)->QVGA(320x240)",176,144, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF(176x144)->VGA (640x480)",176,144, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF (352x288)->SQCIF(128x96)",352,288, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF (352x288)->QCIF(176x144)",352,288, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF (352x288)->QVGA(320x240)",352,288, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF (352x288)->VGA (640x480)",352,288, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA(320x240)->SQCIF(128x96)",320,240, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA(320x240)->QCIF(176x144)",320,240, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA(320x240)->CIF (352x288)",320,240, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA(320x240)->VGA (640x480)",320,240, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA (640x480)->SQCIF(128x96)",640,480, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA (640x480)->QCIF(176x144)",640,480, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA (640x480)->CIF (352x288)",640,480, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA (640x480)->QVGA(320x240)",640,480, 320,240,"bicubic",2);

}

static void mips_resize_bicubic(int isFull, int isVerbose, FILE * fout)
{
  mips_resize_bicubic_gu8(isFull, isVerbose, fout);
  mips_resize_bicubic_gs8(isFull, isVerbose, fout);
  mips_resize_bicubic_gs16(isFull, isVerbose, fout);
}

static void mips_pad(int isFull, int isVerbose, FILE * fout)
{
    // pad
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"padding, SQCIF(128x96)->QCIF(176x144)",128,96, 176,144  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, SQCIF(128x96)->CIF (352x288)",128,96, 352,288  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, SQCIF(128x96)->QVGA(320x240)",128,96, 320,240  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, SQCIF(128x96)->VGA (640x480)" ,128,96, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QCIF(176x144)->CIF (352x288)",176,144, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QCIF(176x144)->QVGA(320x240)",176,144, 320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QCIF(176x144)->VGA (640x480)",176,144, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QVGA(320x240)->CIF (352x288)",320,240, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QVGA(320x240)->VGA (640x480)",320,240, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, CIF (352x288)->VGA (640x480)",352,288, 640,480 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"padding, SQCIF(128x96)->QCIF(176x144)",128,96, 176,144  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, SQCIF(128x96)->CIF (352x288)",128,96, 352,288  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, SQCIF(128x96)->QVGA(320x240)",128,96, 320,240  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, SQCIF(128x96)->VGA (640x480)" ,128,96, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QCIF(176x144)->CIF (352x288)",176,144, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QCIF(176x144)->QVGA(320x240)",176,144, 320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QCIF(176x144)->VGA (640x480)",176,144, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QVGA(320x240)->CIF (352x288)",320,240, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QVGA(320x240)->VGA (640x480)",320,240, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, CIF (352x288)->VGA (640x480)",352,288, 640,480 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"padding, SQCIF(128x96)->QCIF(176x144)",128,96, 176,144 , 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, SQCIF(128x96)->CIF (352x288)",128,96, 352,288 , 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, SQCIF(128x96)->QVGA(320x240)",128,96, 320,240 , 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, SQCIF(128x96)->VGA (640x480)" ,128,96, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QCIF(176x144)->CIF (352x288)",176,144, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QCIF(176x144)->QVGA(320x240)",176,144, 320,240, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QCIF(176x144)->VGA (640x480)",176,144, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QVGA(320x240)->CIF (352x288)",320,240, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QVGA(320x240)->VGA (640x480)",320,240, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, CIF (352x288)->VGA (640x480)",352,288, 640,480, 2);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF(128x96)->QCIF(176x144)",128,96, 176,144 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF(128x96)->CIF (352x288)",128,96, 352,288 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF(128x96)->QVGA(320x240)",128,96, 320,240 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF(128x96)->VGA (640x480)" ,128,96, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QCIF(176x144)->CIF (352x288)",176,144, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QCIF(176x144)->QVGA(320x240)",176,144, 320,240, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QCIF(176x144)->VGA (640x480)",176,144, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QVGA(320x240)->CIF (352x288)",320,240, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QVGA(320x240)->VGA (640x480)",320,240, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, CIF (352x288)->VGA (640x480)",352,288, 640,480, 1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF(128x96)->QCIF(176x144)",128,96, 176,144 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF(128x96)->CIF (352x288)",128,96, 352,288 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF(128x96)->QVGA(320x240)",128,96, 320,240 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF(128x96)->VGA (640x480)" ,128,96, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QCIF(176x144)->CIF (352x288)",176,144, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QCIF(176x144)->QVGA(320x240)",176,144, 320,240, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QCIF(176x144)->VGA (640x480)",176,144, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QVGA(320x240)->CIF (352x288)",320,240, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QVGA(320x240)->VGA (640x480)",320,240, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, CIF (352x288)->VGA (640x480)",352,288, 640,480, 1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF(128x96)->QCIF(176x144)",128,96 , 176,144, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF(128x96)->CIF (352x288)",128,96 , 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF(128x96)->QVGA(320x240)",128,96 , 320,240, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF(128x96)->VGA (640x480)",128,96 , 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QCIF(176x144)->CIF (352x288)",176,144, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QCIF(176x144)->QVGA(320x240)",176,144, 320,240, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QCIF(176x144)->VGA (640x480)",176,144, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QVGA(320x240)->CIF (352x288)",320,240, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QVGA(320x240)->VGA (640x480)",320,240, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, CIF (352x288)->VGA (640x480)",352,288, 640,480, 2);
    // crop
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, QCIF(176x144)->SQCIF(128x96)",176,144,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, CIF (352x288)->SQCIF(128x96)",352,288,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, QVGA(320x240)->SQCIF(128x96)",320,240,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, VGA (640x480)->SQCIF(128x96)",640,480,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, CIF (352x288)->QCIF(176x144)",352,288,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, QVGA(320x240)->QCIF(176x144)",320,240,176,144,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, VGA (640x480)->QCIF(176x144)",640,480,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, CIF (352x288)->QVGA(320x240)",352,288,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, VGA (640x480)->QVGA(320x240)",640,480,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, VGA (640x480)->CIF (352x288)",640,480,352,288,  1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, QCIF(176x144)->SQCIF(128x96)",176,144,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, CIF (352x288)->SQCIF(128x96)",352,288,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, QVGA(320x240)->SQCIF(128x96)",320,240,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, VGA (640x480)->SQCIF(128x96)",640,480,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, CIF (352x288)->QCIF(176x144)",352,288,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, QVGA(320x240)->QCIF(176x144)",320,240,176,144,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, VGA (640x480)->QCIF(176x144)",640,480,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, CIF (352x288)->QVGA(320x240)",352,288,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, VGA (640x480)->QVGA(320x240)",640,480,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, VGA (640x480)->CIF (352x288)",640,480,352,288,  1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"cropping, QCIF(176x144)->SQCIF(128x96)",176,144,128,96 ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"cropping, CIF (352x288)->SQCIF(128x96)",352,288,128,96 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, QVGA(320x240)->SQCIF(128x96)",320,240,128,96 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, VGA (640x480)->SQCIF(128x96)",640,480,128,96 ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"cropping, CIF (352x288)->QCIF(176x144)",352,288,176,144,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, QVGA(320x240)->QCIF(176x144)",320,240,176,144,2);
    PROFILE_IMG_PAD(1     ,isVerbose,imgpad_gs16,0,fout,"cropping, VGA (640x480)->QCIF(176x144)",640,480,176,144,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, CIF (352x288)->QVGA(320x240)",352,288,320,240,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, VGA (640x480)->QVGA(320x240)",640,480,320,240,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, VGA (640x480)->CIF (352x288)",640,480,352,288,2);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"cropping, QCIF(176x144)->SQCIF(128x96)",176,144,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"cropping, CIF (352x288)->SQCIF(128x96)",352,288,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, QVGA(320x240)->SQCIF(128x96)",320,240,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA (640x480)->SQCIF(128x96)",640,480,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"cropping, CIF (352x288)->QCIF(176x144)",352,288,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, QVGA(320x240)->QCIF(176x144)",320,240,176,144 ,1);
    PROFILE_IMG_PAD(1     ,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA (640x480)->QCIF(176x144)",640,480,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, CIF (352x288)->QVGA(320x240)",352,288,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA (640x480)->QVGA(320x240)",640,480,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA (640x480)->CIF (352x288)",640,480,352,288 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"cropping, QCIF(176x144)->SQCIF(128x96)",176,144,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"cropping, CIF (352x288)->SQCIF(128x96)",352,288,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, QVGA(320x240)->SQCIF(128x96)",320,240,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA (640x480)->SQCIF(128x96)",640,480,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"cropping, CIF (352x288)->QCIF(176x144)",352,288,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, QVGA(320x240)->QCIF(176x144)",320,240,176,144 ,1);
    PROFILE_IMG_PAD(1     ,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA (640x480)->QCIF(176x144)",640,480,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, CIF (352x288)->QVGA(320x240)",352,288,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA (640x480)->QVGA(320x240)",640,480,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA (640x480)->CIF (352x288)",640,480,352,288 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"cropping, QCIF(176x144)->SQCIF(128x96)",176,144,128,96  ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"cropping, CIF (352x288)->SQCIF(128x96)",352,288,128,96  ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, QVGA(320x240)->SQCIF(128x96)",320,240,128,96  ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA (640x480)->SQCIF(128x96)",640,480,128,96  ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"cropping, CIF (352x288)->QCIF(176x144)",352,288,176,144 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, QVGA(320x240)->QCIF(176x144)",320,240,176,144 ,2);
    PROFILE_IMG_PAD(1     ,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA (640x480)->QCIF(176x144)",640,480,176,144 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, CIF (352x288)->QVGA(320x240)",352,288,320,240 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA (640x480)->QVGA(320x240)",640,480,320,240 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA (640x480)->CIF (352x288)",640,480,352,288 ,2);
}

void mips_imgresize(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (0==phaseNum || 1==phaseNum) 
  {
    mips_resize_nearest(isFull,isVerbose,fout);
    mips_resize_bilinear(isFull,isVerbose,fout);
    mips_resize_bicubic(isFull,isVerbose,fout);
  }
}

void mips_imgmisc(int phaseNum, int isFull, int isVerbose, FILE * fout)
{
  if (0==phaseNum || 1==phaseNum) 
  {
    mips_hist(isFull,isVerbose,fout);
    mips_norm(isFull,isVerbose,fout);
    mips_interleave(isFull,isVerbose,fout);
    mips_convert(isFull,isVerbose,fout);
    mips_pad(isFull,isVerbose,fout);
  }
} /* mips_imgmisc() */

void mips_imgfft(int phaseNum, int isFull, int isVerbose, FILE * fout)
{ 
  if (0==phaseNum || 1==phaseNum) 
  {
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gu8,fout,"64x64",   64,  64);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"128x128",128, 128);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"256x256",256, 256);
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gu8,fout,"512x512",512, 512);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"SQCIF(128x96)", 128, 96);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"QCIF(176x144)", 176,144);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"QVGA(320x240)", 320,240);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gu8,fout,"CIF (352x288)", 352,288);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gu8,fout,"VGA (640x480)", 640,480);

    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs8,fout,"64x64",   64,  64);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"128x128",128, 128);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"256x256",256, 256);
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs8,fout,"512x512",512, 512);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"SQCIF(128x96)", 128, 96);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"QCIF(176x144)", 176,144);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"QVGA(320x240)", 320,240);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs8,fout,"CIF (352x288)", 352,288);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs8,fout,"VGA (640x480)", 640,480);

    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs16,fout,"64x64",   64,  64);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"128x128",128, 128);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"256x256",256, 256);
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs16,fout,"512x512",512, 512);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"SQCIF(128x96)", 128, 96);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"QCIF(176x144)", 176,144);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"QVGA(320x240)", 320,240);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs16,fout,"CIF (352x288)", 352,288);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs16,fout,"VGA (640x480)", 640,480);

    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gu8,fout,"64x64",   64,  64);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"128x128",128, 128);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"256x256",256, 256);
    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gu8,fout,"512x512",512, 512);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"SQCIF(128x96)", 128, 96);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"QCIF(176x144)", 176,144);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"QVGA(320x240)", 320,240);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gu8,fout,"CIF (352x288)", 352,288);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gu8,fout,"VGA (640x480)", 640,480);

    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs8,fout,"64x64",   64,  64);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"128x128",128, 128);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"256x256",256, 256);
    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs8,fout,"512x512",512, 512);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"SQCIF(128x96)", 128, 96);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"QCIF(176x144)", 176,144);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"QVGA(320x240)", 320,240);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs8,fout,"CIF (352x288)", 352,288);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs8,fout,"VGA (640x480)", 640,480);

    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs16,fout,"64x64",   64,  64);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"128x128",128, 128);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"256x256",256, 256);
    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs16,fout,"512x512",512, 512);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"SQCIF(128x96)", 128, 96);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"QCIF(176x144)", 176,144);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"QVGA(320x240)", 320,240);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs16,fout,"CIF (352x288)", 352,288);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs16,fout,"VGA (640x480)", 640,480);
  }
} /* mips_imgfft() */

