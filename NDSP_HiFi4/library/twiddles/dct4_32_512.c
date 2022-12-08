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

static const complex_fract32 ALIGN(8) dct4_twd_512[]=
{
    {{(int32_t)0X7FFFF621UL,(int32_t)0xFFCDBC0BUL}},
    {{(int32_t)0X7FFFA72CUL,(int32_t)0xFF69343FUL}},
    {{(int32_t)0X7FFF0943UL,(int32_t)0xFF04ACD0UL}},
    {{(int32_t)0X7FFE1C65UL,(int32_t)0xFEA025FDUL}},
    {{(int32_t)0X7FFCE093UL,(int32_t)0xFE3BA002UL}},
    {{(int32_t)0X7FFB55CEUL,(int32_t)0xFDD71B1EUL}},
    {{(int32_t)0X7FF97C18UL,(int32_t)0xFD729790UL}},
    {{(int32_t)0X7FF75370UL,(int32_t)0xFD0E1594UL}},
    {{(int32_t)0X7FF4DBD9UL,(int32_t)0xFCA9956AUL}},
    {{(int32_t)0X7FF21553UL,(int32_t)0xFC45174EUL}},
    {{(int32_t)0X7FEEFFE1UL,(int32_t)0xFBE09B80UL}},
    {{(int32_t)0X7FEB9B85UL,(int32_t)0xFB7C223DUL}},
    {{(int32_t)0X7FE7E841UL,(int32_t)0xFB17ABC2UL}},
    {{(int32_t)0X7FE3E616UL,(int32_t)0xFAB3384FUL}},
    {{(int32_t)0X7FDF9508UL,(int32_t)0xFA4EC821UL}},
    {{(int32_t)0X7FDAF519UL,(int32_t)0xF9EA5B75UL}},
    {{(int32_t)0X7FD6064CUL,(int32_t)0xF985F28AUL}},
    {{(int32_t)0X7FD0C8A3UL,(int32_t)0xF9218D9EUL}},
    {{(int32_t)0X7FCB3C23UL,(int32_t)0xF8BD2CEFUL}},
    {{(int32_t)0X7FC560CFUL,(int32_t)0xF858D0BBUL}},
    {{(int32_t)0X7FBF36AAUL,(int32_t)0xF7F4793EUL}},
    {{(int32_t)0X7FB8BDB8UL,(int32_t)0xF79026B9UL}},
    {{(int32_t)0X7FB1F5FCUL,(int32_t)0xF72BD967UL}},
    {{(int32_t)0X7FAADF7CUL,(int32_t)0xF6C79188UL}},
    {{(int32_t)0X7FA37A3CUL,(int32_t)0xF6634F59UL}},
    {{(int32_t)0X7F9BC640UL,(int32_t)0xF5FF1318UL}},
    {{(int32_t)0X7F93C38CUL,(int32_t)0xF59ADD02UL}},
    {{(int32_t)0X7F8B7227UL,(int32_t)0xF536AD56UL}},
    {{(int32_t)0X7F82D214UL,(int32_t)0xF4D28451UL}},
    {{(int32_t)0X7F79E35AUL,(int32_t)0xF46E6231UL}},
    {{(int32_t)0X7F70A5FEUL,(int32_t)0xF40A4735UL}},
    {{(int32_t)0X7F671A05UL,(int32_t)0xF3A63398UL}},
    {{(int32_t)0X7F5D3F75UL,(int32_t)0xF342279BUL}},
    {{(int32_t)0X7F531655UL,(int32_t)0xF2DE2379UL}},
    {{(int32_t)0X7F489EAAUL,(int32_t)0xF27A2771UL}},
    {{(int32_t)0X7F3DD87CUL,(int32_t)0xF21633C0UL}},
    {{(int32_t)0X7F32C3D1UL,(int32_t)0xF1B248A5UL}},
    {{(int32_t)0X7F2760AFUL,(int32_t)0xF14E665CUL}},
    {{(int32_t)0X7F1BAF1EUL,(int32_t)0xF0EA8D24UL}},
    {{(int32_t)0X7F0FAF25UL,(int32_t)0xF086BD39UL}},
    {{(int32_t)0X7F0360CBUL,(int32_t)0xF022F6DAUL}},
    {{(int32_t)0X7EF6C418UL,(int32_t)0xEFBF3A45UL}},
    {{(int32_t)0X7EE9D914UL,(int32_t)0xEF5B87B5UL}},
    {{(int32_t)0X7EDC9FC6UL,(int32_t)0xEEF7DF6AUL}},
    {{(int32_t)0X7ECF1837UL,(int32_t)0xEE9441A0UL}},
    {{(int32_t)0X7EC14270UL,(int32_t)0xEE30AE96UL}},
    {{(int32_t)0X7EB31E78UL,(int32_t)0xEDCD2687UL}},
    {{(int32_t)0X7EA4AC58UL,(int32_t)0xED69A9B3UL}},
    {{(int32_t)0X7E95EC1AUL,(int32_t)0xED063856UL}},
    {{(int32_t)0X7E86DDC6UL,(int32_t)0xECA2D2ADUL}},
    {{(int32_t)0X7E778166UL,(int32_t)0xEC3F78F6UL}},
    {{(int32_t)0X7E67D703UL,(int32_t)0xEBDC2B6EUL}},
    {{(int32_t)0X7E57DEA7UL,(int32_t)0xEB78EA52UL}},
    {{(int32_t)0X7E47985BUL,(int32_t)0xEB15B5E1UL}},
    {{(int32_t)0X7E37042AUL,(int32_t)0xEAB28E56UL}},
    {{(int32_t)0X7E26221FUL,(int32_t)0xEA4F73EEUL}},
    {{(int32_t)0X7E14F242UL,(int32_t)0xE9EC66E8UL}},
    {{(int32_t)0X7E0374A0UL,(int32_t)0xE9896781UL}},
    {{(int32_t)0X7DF1A942UL,(int32_t)0xE92675F4UL}},
    {{(int32_t)0X7DDF9034UL,(int32_t)0xE8C39280UL}},
    {{(int32_t)0X7DCD2981UL,(int32_t)0xE860BD61UL}},
    {{(int32_t)0X7DBA7534UL,(int32_t)0xE7FDF6D4UL}},
    {{(int32_t)0X7DA77359UL,(int32_t)0xE79B3F16UL}},
    {{(int32_t)0X7D9423FCUL,(int32_t)0xE7389665UL}},
    {{(int32_t)0X7D808728UL,(int32_t)0xE6D5FCFCUL}},
    {{(int32_t)0X7D6C9CE9UL,(int32_t)0xE6737319UL}},
    {{(int32_t)0X7D58654DUL,(int32_t)0xE610F8F9UL}},
    {{(int32_t)0X7D43E05EUL,(int32_t)0xE5AE8ED8UL}},
    {{(int32_t)0X7D2F0E2BUL,(int32_t)0xE54C34F3UL}},
    {{(int32_t)0X7D19EEBFUL,(int32_t)0xE4E9EB87UL}},
    {{(int32_t)0X7D048228UL,(int32_t)0xE487B2D0UL}},
    {{(int32_t)0X7CEEC873UL,(int32_t)0xE4258B0AUL}},
    {{(int32_t)0X7CD8C1AEUL,(int32_t)0xE3C37474UL}},
    {{(int32_t)0X7CC26DE5UL,(int32_t)0xE3616F48UL}},
    {{(int32_t)0X7CABCD28UL,(int32_t)0xE2FF7BC3UL}},
    {{(int32_t)0X7C94DF83UL,(int32_t)0xE29D9A23UL}},
    {{(int32_t)0X7C7DA505UL,(int32_t)0xE23BCAA2UL}},
    {{(int32_t)0X7C661DBCUL,(int32_t)0xE1DA0D7EUL}},
    {{(int32_t)0X7C4E49B7UL,(int32_t)0xE17862F3UL}},
    {{(int32_t)0X7C362904UL,(int32_t)0xE116CB3DUL}},
    {{(int32_t)0X7C1DBBB3UL,(int32_t)0xE0B54698UL}},
    {{(int32_t)0X7C0501D2UL,(int32_t)0xE053D541UL}},
    {{(int32_t)0X7BEBFB70UL,(int32_t)0xDFF27773UL}},
    {{(int32_t)0X7BD2A89EUL,(int32_t)0xDF912D6BUL}},
    {{(int32_t)0X7BB9096BUL,(int32_t)0xDF2FF764UL}},
    {{(int32_t)0X7B9F1DE6UL,(int32_t)0xDECED59BUL}},
    {{(int32_t)0X7B84E61FUL,(int32_t)0xDE6DC84BUL}},
    {{(int32_t)0X7B6A6227UL,(int32_t)0xDE0CCFB1UL}},
    {{(int32_t)0X7B4F920EUL,(int32_t)0xDDABEC08UL}},
    {{(int32_t)0X7B3475E5UL,(int32_t)0xDD4B1D8CUL}},
    {{(int32_t)0X7B190DBCUL,(int32_t)0xDCEA6478UL}},
    {{(int32_t)0X7AFD59A4UL,(int32_t)0xDC89C109UL}},
    {{(int32_t)0X7AE159AEUL,(int32_t)0xDC293379UL}},
    {{(int32_t)0X7AC50DECUL,(int32_t)0xDBC8BC06UL}},
    {{(int32_t)0X7AA8766FUL,(int32_t)0xDB685AE9UL}},
    {{(int32_t)0X7A8B9348UL,(int32_t)0xDB08105EUL}},
    {{(int32_t)0X7A6E648AUL,(int32_t)0xDAA7DCA1UL}},
    {{(int32_t)0X7A50EA47UL,(int32_t)0xDA47BFEEUL}},
    {{(int32_t)0X7A332490UL,(int32_t)0xD9E7BA7FUL}},
    {{(int32_t)0X7A151378UL,(int32_t)0xD987CC90UL}},
    {{(int32_t)0X79F6B711UL,(int32_t)0xD927F65BUL}},
    {{(int32_t)0X79D80F6FUL,(int32_t)0xD8C8381DUL}},
    {{(int32_t)0X79B91CA4UL,(int32_t)0xD868920FUL}},
    {{(int32_t)0X7999DEC4UL,(int32_t)0xD809046EUL}},
    {{(int32_t)0X797A55E0UL,(int32_t)0xD7A98F73UL}},
    {{(int32_t)0X795A820EUL,(int32_t)0xD74A335BUL}},
    {{(int32_t)0X793A6361UL,(int32_t)0xD6EAF05FUL}},
    {{(int32_t)0X7919F9ECUL,(int32_t)0xD68BC6BAUL}},
    {{(int32_t)0X78F945C3UL,(int32_t)0xD62CB6A8UL}},
    {{(int32_t)0X78D846FBUL,(int32_t)0xD5CDC062UL}},
    {{(int32_t)0X78B6FDA8UL,(int32_t)0xD56EE424UL}},
    {{(int32_t)0X789569DFUL,(int32_t)0xD5102228UL}},
    {{(int32_t)0X78738BB3UL,(int32_t)0xD4B17AA8UL}},
    {{(int32_t)0X7851633BUL,(int32_t)0xD452EDDFUL}},
    {{(int32_t)0X782EF08BUL,(int32_t)0xD3F47C06UL}},
    {{(int32_t)0X780C33B8UL,(int32_t)0xD396255AUL}},
    {{(int32_t)0X77E92CD9UL,(int32_t)0xD337EA12UL}},
    {{(int32_t)0X77C5DC01UL,(int32_t)0xD2D9CA6AUL}},
    {{(int32_t)0X77A24148UL,(int32_t)0xD27BC69CUL}},
    {{(int32_t)0X777E5CC3UL,(int32_t)0xD21DDEE2UL}},
    {{(int32_t)0X775A2E89UL,(int32_t)0xD1C01375UL}},
    {{(int32_t)0X7735B6AFUL,(int32_t)0xD1626490UL}},
    {{(int32_t)0X7710F54CUL,(int32_t)0xD104D26BUL}},
    {{(int32_t)0X76EBEA77UL,(int32_t)0xD0A75D42UL}},
    {{(int32_t)0X76C69647UL,(int32_t)0xD04A054EUL}},
    {{(int32_t)0X76A0F8D2UL,(int32_t)0xCFECCAC7UL}},
    {{(int32_t)0X767B1231UL,(int32_t)0xCF8FADE9UL}},
    {{(int32_t)0X7654E279UL,(int32_t)0xCF32AEEBUL}},
    {{(int32_t)0X5A5EE79AUL,(int32_t)0xA55A025BUL}},
    {{(int32_t)0X5A1799D1UL,(int32_t)0xA513243BUL}},
    {{(int32_t)0X59D01475UL,(int32_t)0xA4CC7E32UL}},
    {{(int32_t)0X598857B2UL,(int32_t)0xA486106AUL}},
    {{(int32_t)0X594063B5UL,(int32_t)0xA43FDB10UL}},
    {{(int32_t)0X58F838A9UL,(int32_t)0xA3F9DE4EUL}},
    {{(int32_t)0X58AFD6BDUL,(int32_t)0xA3B41A50UL}},
    {{(int32_t)0X58673E1BUL,(int32_t)0xA36E8F41UL}},
    {{(int32_t)0X581E6EF1UL,(int32_t)0xA3293D4BUL}},
    {{(int32_t)0X57D5696DUL,(int32_t)0xA2E4249BUL}},
    {{(int32_t)0X578C2DBAUL,(int32_t)0xA29F4559UL}},
    {{(int32_t)0X5742BC06UL,(int32_t)0xA25A9FB1UL}},
    {{(int32_t)0X56F9147EUL,(int32_t)0xA21633CDUL}},
    {{(int32_t)0X56AF3750UL,(int32_t)0xA1D201D7UL}},
    {{(int32_t)0X566524AAUL,(int32_t)0xA18E09FAUL}},
    {{(int32_t)0X561ADCB9UL,(int32_t)0xA14A4C5EUL}},
    {{(int32_t)0X55D05FAAUL,(int32_t)0xA106C92FUL}},
    {{(int32_t)0X5585ADADUL,(int32_t)0xA0C38095UL}},
    {{(int32_t)0X553AC6EEUL,(int32_t)0xA08072BAUL}},
    {{(int32_t)0X54EFAB9CUL,(int32_t)0xA03D9FC8UL}},
    {{(int32_t)0X54A45BE6UL,(int32_t)0x9FFB07E7UL}},
    {{(int32_t)0X5458D7F9UL,(int32_t)0x9FB8AB41UL}},
    {{(int32_t)0X540D2005UL,(int32_t)0x9F7689FFUL}},
    {{(int32_t)0X53C13439UL,(int32_t)0x9F34A449UL}},
    {{(int32_t)0X537514C2UL,(int32_t)0x9EF2FA49UL}},
    {{(int32_t)0X5328C1D0UL,(int32_t)0x9EB18C26UL}},
    {{(int32_t)0X52DC3B92UL,(int32_t)0x9E705A09UL}},
    {{(int32_t)0X528F8238UL,(int32_t)0x9E2F641BUL}},
    {{(int32_t)0X524295F0UL,(int32_t)0x9DEEAA82UL}},
    {{(int32_t)0X51F576EAUL,(int32_t)0x9DAE2D68UL}},
    {{(int32_t)0X51A82555UL,(int32_t)0x9D6DECF4UL}},
    {{(int32_t)0X515AA162UL,(int32_t)0x9D2DE94DUL}},
    {{(int32_t)0X510CEB40UL,(int32_t)0x9CEE229CUL}},
    {{(int32_t)0X50BF031FUL,(int32_t)0x9CAE9907UL}},
    {{(int32_t)0X5070E92FUL,(int32_t)0x9C6F4CB6UL}},
    {{(int32_t)0X50229DA1UL,(int32_t)0x9C303DCFUL}},
    {{(int32_t)0X4FD420A4UL,(int32_t)0x9BF16C7AUL}},
    {{(int32_t)0X4F857269UL,(int32_t)0x9BB2D8DEUL}},
    {{(int32_t)0X4F369320UL,(int32_t)0x9B748320UL}},
    {{(int32_t)0X4EE782FBUL,(int32_t)0x9B366B68UL}},
    {{(int32_t)0X4E984229UL,(int32_t)0x9AF891DBUL}},
    {{(int32_t)0X4E48D0DDUL,(int32_t)0x9ABAF6A1UL}},
    {{(int32_t)0X4DF92F46UL,(int32_t)0x9A7D99DEUL}},
    {{(int32_t)0X4DA95D96UL,(int32_t)0x9A407BB9UL}},
    {{(int32_t)0X4D595BFEUL,(int32_t)0x9A039C57UL}},
    {{(int32_t)0X4D092AB0UL,(int32_t)0x99C6FBDEUL}},
    {{(int32_t)0X4CB8C9DDUL,(int32_t)0x998A9A74UL}},
    {{(int32_t)0X4C6839B7UL,(int32_t)0x994E783DUL}},
    {{(int32_t)0X4C177A6EUL,(int32_t)0x9912955FUL}},
    {{(int32_t)0X4BC68C36UL,(int32_t)0x98D6F1FEUL}},
    {{(int32_t)0X4B756F40UL,(int32_t)0x989B8E40UL}},
    {{(int32_t)0X4B2423BEUL,(int32_t)0x98606A49UL}},
    {{(int32_t)0X4AD2A9E2UL,(int32_t)0x9825863DUL}},
    {{(int32_t)0X4A8101DEUL,(int32_t)0x97EAE242UL}},
    {{(int32_t)0X4A2F2BE6UL,(int32_t)0x97B07E7AUL}},
    {{(int32_t)0X49DD282AUL,(int32_t)0x97765B0AUL}},
    {{(int32_t)0X498AF6DFUL,(int32_t)0x973C7817UL}},
    {{(int32_t)0X49389836UL,(int32_t)0x9702D5C3UL}},
    {{(int32_t)0X48E60C62UL,(int32_t)0x96C97432UL}},
    {{(int32_t)0X48935397UL,(int32_t)0x96905388UL}},
    {{(int32_t)0X48406E08UL,(int32_t)0x965773E7UL}},
    {{(int32_t)0X47ED5BE6UL,(int32_t)0x961ED574UL}},
    {{(int32_t)0X479A1D67UL,(int32_t)0x95E67850UL}},
    {{(int32_t)0X4746B2BCUL,(int32_t)0x95AE5C9FUL}},
    {{(int32_t)0X46F31C1AUL,(int32_t)0x95768283UL}},
    {{(int32_t)0X469F59B4UL,(int32_t)0x953EEA1EUL}},
    {{(int32_t)0X464B6BBEUL,(int32_t)0x95079394UL}},
    {{(int32_t)0X45F7526BUL,(int32_t)0x94D07F05UL}},
    {{(int32_t)0X45A30DF0UL,(int32_t)0x9499AC95UL}},
    {{(int32_t)0X454E9E80UL,(int32_t)0x94631C65UL}},
    {{(int32_t)0X44FA0450UL,(int32_t)0x942CCE96UL}},
    {{(int32_t)0X44A53F93UL,(int32_t)0x93F6C34AUL}},
    {{(int32_t)0X4450507EUL,(int32_t)0x93C0FAA3UL}},
    {{(int32_t)0X43FB3746UL,(int32_t)0x938B74C1UL}},
    {{(int32_t)0X43A5F41EUL,(int32_t)0x935631C5UL}},
    {{(int32_t)0X4350873CUL,(int32_t)0x932131D1UL}},
    {{(int32_t)0X42FAF0D4UL,(int32_t)0x92EC7505UL}},
    {{(int32_t)0X42A5311BUL,(int32_t)0x92B7FB82UL}},
    {{(int32_t)0X424F4845UL,(int32_t)0x9283C568UL}},
    {{(int32_t)0X41F93689UL,(int32_t)0x924FD2D7UL}},
    {{(int32_t)0X41A2FC1AUL,(int32_t)0x921C23EFUL}},
    {{(int32_t)0X414C992FUL,(int32_t)0x91E8B8D0UL}},
    {{(int32_t)0X40F60DFBUL,(int32_t)0x91B5919AUL}},
    {{(int32_t)0X409F5AB6UL,(int32_t)0x9182AE6DUL}},
    {{(int32_t)0X40487F94UL,(int32_t)0x91500F67UL}},
    {{(int32_t)0X3FF17CCAUL,(int32_t)0x911DB4A9UL}},
    {{(int32_t)0X3F9A5290UL,(int32_t)0x90EB9E50UL}},
    {{(int32_t)0X3F430119UL,(int32_t)0x90B9CC7DUL}},
    {{(int32_t)0X3EEB889CUL,(int32_t)0x90883F4DUL}},
    {{(int32_t)0X3E93E950UL,(int32_t)0x9056F6DFUL}},
    {{(int32_t)0X3E3C2369UL,(int32_t)0x9025F352UL}},
    {{(int32_t)0X3DE4371FUL,(int32_t)0x8FF534C4UL}},
    {{(int32_t)0X3D8C24A8UL,(int32_t)0x8FC4BB53UL}},
    {{(int32_t)0X3D33EC39UL,(int32_t)0x8F94871DUL}},
    {{(int32_t)0X3CDB8E09UL,(int32_t)0x8F649840UL}},
    {{(int32_t)0X3C830A50UL,(int32_t)0x8F34EED8UL}},
    {{(int32_t)0X3C2A6142UL,(int32_t)0x8F058B04UL}},
    {{(int32_t)0X3BD19318UL,(int32_t)0x8ED66CE1UL}},
    {{(int32_t)0X3B78A007UL,(int32_t)0x8EA7948CUL}},
    {{(int32_t)0X3B1F8848UL,(int32_t)0x8E790222UL}},
    {{(int32_t)0X3AC64C0FUL,(int32_t)0x8E4AB5BFUL}},
    {{(int32_t)0X3A6CEB96UL,(int32_t)0x8E1CAF80UL}},
    {{(int32_t)0X3A136712UL,(int32_t)0x8DEEEF82UL}},
    {{(int32_t)0X39B9BEBCUL,(int32_t)0x8DC175E0UL}},
    {{(int32_t)0X395FF2C9UL,(int32_t)0x8D9442B8UL}},
    {{(int32_t)0X39060373UL,(int32_t)0x8D675623UL}},
    {{(int32_t)0X38ABF0EFUL,(int32_t)0x8D3AB03FUL}},
    {{(int32_t)0X3851BB77UL,(int32_t)0x8D0E5127UL}},
    {{(int32_t)0X37F76341UL,(int32_t)0x8CE238F6UL}},
    {{(int32_t)0X379CE885UL,(int32_t)0x8CB667C8UL}},
    {{(int32_t)0X37424B7BUL,(int32_t)0x8C8ADDB7UL}},
    {{(int32_t)0X36E78C5BUL,(int32_t)0x8C5F9ADEUL}},
    {{(int32_t)0X368CAB5CUL,(int32_t)0x8C349F58UL}},
    {{(int32_t)0X3631A8B8UL,(int32_t)0x8C09EB40UL}},
    {{(int32_t)0X35D684A6UL,(int32_t)0x8BDF7EB0UL}},
    {{(int32_t)0X357B3F5DUL,(int32_t)0x8BB559C1UL}},
    {{(int32_t)0X351FD918UL,(int32_t)0x8B8B7C8FUL}},
    {{(int32_t)0X34C4520DUL,(int32_t)0x8B61E733UL}},
    {{(int32_t)0X3468AA76UL,(int32_t)0x8B3899C6UL}},
    {{(int32_t)0X340CE28BUL,(int32_t)0x8B0F9462UL}},
    {{(int32_t)0X33B0FA84UL,(int32_t)0x8AE6D720UL}},
    {{(int32_t)0X3354F29BUL,(int32_t)0x8ABE6219UL}},
    {{(int32_t)0X32F8CB07UL,(int32_t)0x8A963567UL}},
    {{(int32_t)0X329C8402UL,(int32_t)0x8A6E5123UL}},
    {{(int32_t)0X32401DC6UL,(int32_t)0x8A46B564UL}},
    {{(int32_t)0X31E39889UL,(int32_t)0x8A1F6243UL}},
    {{(int32_t)0X3186F487UL,(int32_t)0x89F857D8UL}},
    {{(int32_t)0X312A31F8UL,(int32_t)0x89D1963CUL}},
};

static const complex_fract32 dct3_512[]=
{
    {{(int32_t)0X7FFFFFFF,(int32_t)0000000000UL}},
    {{(int32_t)0X7FFF6216UL,(int32_t)0x00C90F88UL}},
    {{(int32_t)0X7FFD885AUL,(int32_t)0x01921D20UL}},
    {{(int32_t)0X7FFA72D1UL,(int32_t)0x025B26D7UL}},
    {{(int32_t)0X7FF62182UL,(int32_t)0x03242ABFUL}},
    {{(int32_t)0X7FF09478UL,(int32_t)0x03ED26E6UL}},
    {{(int32_t)0X7FE9CBC0UL,(int32_t)0x04B6195DUL}},
    {{(int32_t)0X7FE1C76BUL,(int32_t)0x057F0035UL}},
    {{(int32_t)0X7FD8878EUL,(int32_t)0x0647D97CUL}},
    {{(int32_t)0X7FCE0C3EUL,(int32_t)0x0710A345UL}},
    {{(int32_t)0X7FC25596UL,(int32_t)0x07D95B9EUL}},
    {{(int32_t)0X7FB563B3UL,(int32_t)0x08A2009AUL}},
    {{(int32_t)0X7FA736B4UL,(int32_t)0x096A9049UL}},
    {{(int32_t)0X7F97CEBDUL,(int32_t)0x0A3308BDUL}},
    {{(int32_t)0X7F872BF3UL,(int32_t)0x0AFB6805UL}},
    {{(int32_t)0X7F754E80UL,(int32_t)0x0BC3AC35UL}},
    {{(int32_t)0X7F62368FUL,(int32_t)0x0C8BD35EUL}},
    {{(int32_t)0X7F4DE451UL,(int32_t)0x0D53DB92UL}},
    {{(int32_t)0X7F3857F6UL,(int32_t)0x0E1BC2E4UL}},
    {{(int32_t)0X7F2191B4UL,(int32_t)0x0EE38766UL}},
    {{(int32_t)0X7F0991C4UL,(int32_t)0x0FAB272BUL}},
    {{(int32_t)0X7EF05860UL,(int32_t)0x1072A048UL}},
    {{(int32_t)0X7ED5E5C6UL,(int32_t)0x1139F0CFUL}},
    {{(int32_t)0X7EBA3A39UL,(int32_t)0x120116D5UL}},
    {{(int32_t)0X7E9D55FCUL,(int32_t)0x12C8106FUL}},
    {{(int32_t)0X7E7F3957UL,(int32_t)0x138EDBB1UL}},
    {{(int32_t)0X7E5FE493UL,(int32_t)0x145576B1UL}},
    {{(int32_t)0X7E3F57FFUL,(int32_t)0x151BDF86UL}},
    {{(int32_t)0X7E1D93EAUL,(int32_t)0x15E21445UL}},
    {{(int32_t)0X7DFA98A8UL,(int32_t)0x16A81305UL}},
    {{(int32_t)0X7DD6668FUL,(int32_t)0x176DD9DEUL}},
    {{(int32_t)0X7DB0FDF8UL,(int32_t)0x183366E9UL}},
    {{(int32_t)0X7D8A5F40UL,(int32_t)0x18F8B83CUL}},
    {{(int32_t)0X7D628AC6UL,(int32_t)0x19BDCBF3UL}},
    {{(int32_t)0X7D3980ECUL,(int32_t)0x1A82A026UL}},
    {{(int32_t)0X7D0F4218UL,(int32_t)0x1B4732EFUL}},
    {{(int32_t)0X7CE3CEB2UL,(int32_t)0x1C0B826AUL}},
    {{(int32_t)0X7CB72724UL,(int32_t)0x1CCF8CB3UL}},
    {{(int32_t)0X7C894BDEUL,(int32_t)0x1D934FE5UL}},
    {{(int32_t)0X7C5A3D50UL,(int32_t)0x1E56CA1EUL}},
    {{(int32_t)0X7C29FBEEUL,(int32_t)0x1F19F97BUL}},
    {{(int32_t)0X7BF88830UL,(int32_t)0x1FDCDC1BUL}},
    {{(int32_t)0X7BC5E290UL,(int32_t)0x209F701CUL}},
    {{(int32_t)0X7B920B89UL,(int32_t)0x2161B3A0UL}},
    {{(int32_t)0X7B5D039EUL,(int32_t)0x2223A4C5UL}},
    {{(int32_t)0X7B26CB4FUL,(int32_t)0x22E541AFUL}},
    {{(int32_t)0X7AEF6323UL,(int32_t)0x23A6887FUL}},
    {{(int32_t)0X7AB6CBA4UL,(int32_t)0x24677758UL}},
    {{(int32_t)0X7A7D055BUL,(int32_t)0x25280C5EUL}},
    {{(int32_t)0X7A4210D8UL,(int32_t)0x25E845B6UL}},
    {{(int32_t)0X7A05EEADUL,(int32_t)0x26A82186UL}},
    {{(int32_t)0X79C89F6EUL,(int32_t)0x27679DF4UL}},
    {{(int32_t)0X798A23B1UL,(int32_t)0x2826B928UL}},
    {{(int32_t)0X794A7C12UL,(int32_t)0x28E5714BUL}},
    {{(int32_t)0X7909A92DUL,(int32_t)0x29A3C485UL}},
    {{(int32_t)0X78C7ABA2UL,(int32_t)0x2A61B101UL}},
    {{(int32_t)0X78848414UL,(int32_t)0x2B1F34EBUL}},
    {{(int32_t)0X78403329UL,(int32_t)0x2BDC4E6FUL}},
    {{(int32_t)0X77FAB989UL,(int32_t)0x2C98FBBAUL}},
    {{(int32_t)0X77B417DFUL,(int32_t)0x2D553AFCUL}},
    {{(int32_t)0X776C4EDBUL,(int32_t)0x2E110A62UL}},
    {{(int32_t)0X77235F2DUL,(int32_t)0x2ECC681EUL}},
    {{(int32_t)0X76D94989UL,(int32_t)0x2F875262UL}},
    {{(int32_t)0X768E0EA6UL,(int32_t)0x3041C761UL}},
    {{(int32_t)0X7641AF3DUL,(int32_t)0x30FBC54DUL}},
    {{(int32_t)0X75F42C0BUL,(int32_t)0x31B54A5EUL}},
    {{(int32_t)0X75A585CFUL,(int32_t)0x326E54C7UL}},
    {{(int32_t)0X7555BD4CUL,(int32_t)0x3326E2C3UL}},
    {{(int32_t)0X7504D345UL,(int32_t)0x33DEF287UL}},
    {{(int32_t)0X74B2C884UL,(int32_t)0x34968250UL}},
    {{(int32_t)0X745F9DD1UL,(int32_t)0x354D9057UL}},
    {{(int32_t)0X740B53FBUL,(int32_t)0x36041AD9UL}},
    {{(int32_t)0X73B5EBD1UL,(int32_t)0x36BA2014UL}},
    {{(int32_t)0X735F6626UL,(int32_t)0x376F9E46UL}},
    {{(int32_t)0X7307C3D0UL,(int32_t)0x382493B0UL}},
    {{(int32_t)0X72AF05A7UL,(int32_t)0x38D8FE93UL}},
    {{(int32_t)0X72552C85UL,(int32_t)0x398CDD32UL}},
    {{(int32_t)0X71FA3949UL,(int32_t)0x3A402DD2UL}},
    {{(int32_t)0X719E2CD2UL,(int32_t)0x3AF2EEB7UL}},
    {{(int32_t)0X71410805UL,(int32_t)0x3BA51E29UL}},
    {{(int32_t)0X70E2CBC6UL,(int32_t)0x3C56BA70UL}},
    {{(int32_t)0X708378FFUL,(int32_t)0x3D07C1D6UL}},
    {{(int32_t)0X7023109AUL,(int32_t)0x3DB832A6UL}},
    {{(int32_t)0X6FC19385UL,(int32_t)0x3E680B2CUL}},
    {{(int32_t)0X6F5F02B2UL,(int32_t)0x3F1749B8UL}},
    {{(int32_t)0X6EFB5F12UL,(int32_t)0x3FC5EC98UL}},
    {{(int32_t)0X6E96A99DUL,(int32_t)0x4073F21DUL}},
    {{(int32_t)0X6E30E34AUL,(int32_t)0x4121589BUL}},
    {{(int32_t)0X6DCA0D14UL,(int32_t)0x41CE1E65UL}},
    {{(int32_t)0X6D6227FAUL,(int32_t)0x427A41D0UL}},
    {{(int32_t)0X6CF934FCUL,(int32_t)0x4325C135UL}},
    {{(int32_t)0X6C8F351CUL,(int32_t)0x43D09AEDUL}},
    {{(int32_t)0X6C242960UL,(int32_t)0x447ACD50UL}},
    {{(int32_t)0X6BB812D1UL,(int32_t)0x452456BDUL}},
    {{(int32_t)0X6B4AF279UL,(int32_t)0x45CD358FUL}},
    {{(int32_t)0X6ADCC964UL,(int32_t)0x46756828UL}},
    {{(int32_t)0X6A6D98A4UL,(int32_t)0x471CECE7UL}},
    {{(int32_t)0X69FD614AUL,(int32_t)0x47C3C22FUL}},
    {{(int32_t)0X698C246CUL,(int32_t)0x4869E665UL}},
    {{(int32_t)0X6919E320UL,(int32_t)0x490F57EEUL}},
    {{(int32_t)0X68A69E81UL,(int32_t)0x49B41533UL}},
    {{(int32_t)0X683257ABUL,(int32_t)0x4A581C9EUL}},
    {{(int32_t)0X67BD0FBDUL,(int32_t)0x4AFB6C98UL}},
    {{(int32_t)0X6746C7D8UL,(int32_t)0x4B9E0390UL}},
    {{(int32_t)0X66CF8120UL,(int32_t)0x4C3FDFF4UL}},
    {{(int32_t)0X66573CBBUL,(int32_t)0x4CE10034UL}},
    {{(int32_t)0X65DDFBD3UL,(int32_t)0x4D8162C4UL}},
    {{(int32_t)0X6563BF92UL,(int32_t)0x4E210617UL}},
    {{(int32_t)0X64E88926UL,(int32_t)0x4EBFE8A5UL}},
    {{(int32_t)0X646C59BFUL,(int32_t)0x4F5E08E3UL}},
    {{(int32_t)0X63EF3290UL,(int32_t)0x4FFB654DUL}},
    {{(int32_t)0X637114CCUL,(int32_t)0x5097FC5EUL}},
    {{(int32_t)0X62F201ACUL,(int32_t)0x5133CC94UL}},
    {{(int32_t)0X6271FA69UL,(int32_t)0x51CED46EUL}},
    {{(int32_t)0X61F1003FUL,(int32_t)0x5269126EUL}},
    {{(int32_t)0X616F146CUL,(int32_t)0x53028518UL}},
    {{(int32_t)0X60EC3830UL,(int32_t)0x539B2AF0UL}},
    {{(int32_t)0X60686CCFUL,(int32_t)0x5433027DUL}},
    {{(int32_t)0X5FE3B38DUL,(int32_t)0x54CA0A4BUL}},
    {{(int32_t)0X5F5E0DB3UL,(int32_t)0x556040E2UL}},
    {{(int32_t)0X5ED77C8AUL,(int32_t)0x55F5A4D2UL}},
    {{(int32_t)0X5E50015DUL,(int32_t)0x568A34A9UL}},
    {{(int32_t)0X5DC79D7CUL,(int32_t)0x571DEEFAUL}},
    {{(int32_t)0X5D3E5237UL,(int32_t)0x57B0D256UL}},
    {{(int32_t)0X5CB420E0UL,(int32_t)0x5842DD54UL}},
    {{(int32_t)0X5C290ACCUL,(int32_t)0x58D40E8CUL}},
    {{(int32_t)0X5B9D1154UL,(int32_t)0x59646498UL}},
    {{(int32_t)0X5B1035CFUL,(int32_t)0x59F3DE12UL}},
};

static const complex_fract32 rfft_256[]=
{
    {{(int32_t)0X7FF62182UL,(int32_t)0x03242ABFUL}},
    {{(int32_t)0X7FD8878EUL,(int32_t)0x0647D97CUL}},
    {{(int32_t)0X7FA736B4UL,(int32_t)0x096A9049UL}},
    {{(int32_t)0X7F62368FUL,(int32_t)0x0C8BD35EUL}},
    {{(int32_t)0X7F0991C4UL,(int32_t)0x0FAB272BUL}},
    {{(int32_t)0X7E9D55FCUL,(int32_t)0x12C8106FUL}},
    {{(int32_t)0X7E1D93EAUL,(int32_t)0x15E21445UL}},
    {{(int32_t)0X7D8A5F40UL,(int32_t)0x18F8B83CUL}},
    {{(int32_t)0X7CE3CEB2UL,(int32_t)0x1C0B826AUL}},
    {{(int32_t)0X7C29FBEEUL,(int32_t)0x1F19F97BUL}},
    {{(int32_t)0X7B5D039EUL,(int32_t)0x2223A4C5UL}},
    {{(int32_t)0X7A7D055BUL,(int32_t)0x25280C5EUL}},
    {{(int32_t)0X798A23B1UL,(int32_t)0x2826B928UL}},
    {{(int32_t)0X78848414UL,(int32_t)0x2B1F34EBUL}},
    {{(int32_t)0X776C4EDBUL,(int32_t)0x2E110A62UL}},
    {{(int32_t)0X7641AF3DUL,(int32_t)0x30FBC54DUL}},
    {{(int32_t)0X7504D345UL,(int32_t)0x33DEF287UL}},
    {{(int32_t)0X73B5EBD1UL,(int32_t)0x36BA2014UL}},
    {{(int32_t)0X72552C85UL,(int32_t)0x398CDD32UL}},
    {{(int32_t)0X70E2CBC6UL,(int32_t)0x3C56BA70UL}},
    {{(int32_t)0X6F5F02B2UL,(int32_t)0x3F1749B8UL}},
    {{(int32_t)0X6DCA0D14UL,(int32_t)0x41CE1E65UL}},
    {{(int32_t)0X6C242960UL,(int32_t)0x447ACD50UL}},
    {{(int32_t)0X6A6D98A4UL,(int32_t)0x471CECE7UL}},
    {{(int32_t)0X68A69E81UL,(int32_t)0x49B41533UL}},
    {{(int32_t)0X66CF8120UL,(int32_t)0x4C3FDFF4UL}},
    {{(int32_t)0X64E88926UL,(int32_t)0x4EBFE8A5UL}},
    {{(int32_t)0X62F201ACUL,(int32_t)0x5133CC94UL}},
    {{(int32_t)0X60EC3830UL,(int32_t)0x539B2AF0UL}},
    {{(int32_t)0X5ED77C8AUL,(int32_t)0x55F5A4D2UL}},
    {{(int32_t)0X5CB420E0UL,(int32_t)0x5842DD54UL}},
    {{(int32_t)0X5A82799AUL,(int32_t)0x5A82799AUL}},
    {{(int32_t)0X5842DD54UL,(int32_t)0x5CB420E0UL}},
    {{(int32_t)0X55F5A4D2UL,(int32_t)0x5ED77C8AUL}},
    {{(int32_t)0X539B2AF0UL,(int32_t)0x60EC3830UL}},
    {{(int32_t)0X5133CC94UL,(int32_t)0x62F201ACUL}},
    {{(int32_t)0X4EBFE8A5UL,(int32_t)0x64E88926UL}},
    {{(int32_t)0X4C3FDFF4UL,(int32_t)0x66CF8120UL}},
    {{(int32_t)0X49B41533UL,(int32_t)0x68A69E81UL}},
    {{(int32_t)0X471CECE7UL,(int32_t)0x6A6D98A4UL}},
    {{(int32_t)0X447ACD50UL,(int32_t)0x6C242960UL}},
    {{(int32_t)0X41CE1E65UL,(int32_t)0x6DCA0D14UL}},
    {{(int32_t)0X3F1749B8UL,(int32_t)0x6F5F02B2UL}},
    {{(int32_t)0X3C56BA70UL,(int32_t)0x70E2CBC6UL}},
    {{(int32_t)0X398CDD32UL,(int32_t)0x72552C85UL}},
    {{(int32_t)0X36BA2014UL,(int32_t)0x73B5EBD1UL}},
    {{(int32_t)0X33DEF287UL,(int32_t)0x7504D345UL}},
    {{(int32_t)0X30FBC54DUL,(int32_t)0x7641AF3DUL}},
    {{(int32_t)0X2E110A62UL,(int32_t)0x776C4EDBUL}},
    {{(int32_t)0X2B1F34EBUL,(int32_t)0x78848414UL}},
    {{(int32_t)0X2826B928UL,(int32_t)0x798A23B1UL}},
    {{(int32_t)0X25280C5EUL,(int32_t)0x7A7D055BUL}},
    {{(int32_t)0X2223A4C5UL,(int32_t)0x7B5D039EUL}},
    {{(int32_t)0X1F19F97BUL,(int32_t)0x7C29FBEEUL}},
    {{(int32_t)0X1C0B826AUL,(int32_t)0x7CE3CEB2UL}},
    {{(int32_t)0X18F8B83CUL,(int32_t)0x7D8A5F40UL}},
    {{(int32_t)0X15E21445UL,(int32_t)0x7E1D93EAUL}},
    {{(int32_t)0X12C8106FUL,(int32_t)0x7E9D55FCUL}},
    {{(int32_t)0X0FAB272BUL,(int32_t)0x7F0991C4UL}},
    {{(int32_t)0X0C8BD35EUL,(int32_t)0x7F62368FUL}},
    {{(int32_t)0X096A9049UL,(int32_t)0x7FA736B4UL}},
    {{(int32_t)0X0647D97CUL,(int32_t)0x7FD8878EUL}},
    {{(int32_t)0X03242ABFUL,(int32_t)0x7FF62182UL}},
};

static const complex_fract32 fft_128[]=
{
    {{(int32_t)0X7FFFFFFF,(int32_t)0000000000UL}},
    {{(int32_t)0X7FFFFFFF,(int32_t)0000000000UL}},
    {{(int32_t)0X7FFFFFFF,(int32_t)0000000000UL}},
    {{(int32_t)0X7FD8878EUL,(int32_t)0xF9B82684UL}},
    {{(int32_t)0X7F62368FUL,(int32_t)0xF3742CA2UL}},
    {{(int32_t)0X7E9D55FCUL,(int32_t)0xED37EF91UL}},
    {{(int32_t)0X7F62368FUL,(int32_t)0xF3742CA2UL}},
    {{(int32_t)0X7D8A5F40UL,(int32_t)0xE70747C4UL}},
    {{(int32_t)0X7A7D055BUL,(int32_t)0xDAD7F3A2UL}},
    {{(int32_t)0X7E9D55FCUL,(int32_t)0xED37EF91UL}},
    {{(int32_t)0X7A7D055BUL,(int32_t)0xDAD7F3A2UL}},
    {{(int32_t)0X73B5EBD1UL,(int32_t)0xC945DFECUL}},
    {{(int32_t)0X7D8A5F40UL,(int32_t)0xE70747C4UL}},
    {{(int32_t)0X7641AF3DUL,(int32_t)0xCF043AB3UL}},
    {{(int32_t)0X6A6D98A4UL,(int32_t)0xB8E31319UL}},
    {{(int32_t)0X7C29FBEEUL,(int32_t)0xE0E60685UL}},
    {{(int32_t)0X70E2CBC6UL,(int32_t)0xC3A94590UL}},
    {{(int32_t)0X5ED77C8AUL,(int32_t)0xAA0A5B2EUL}},
    {{(int32_t)0X7A7D055BUL,(int32_t)0xDAD7F3A2UL}},
    {{(int32_t)0X6A6D98A4UL,(int32_t)0xB8E31319UL}},
    {{(int32_t)0X5133CC94UL,(int32_t)0x9D0DFE54UL}},
    {{(int32_t)0X78848414UL,(int32_t)0xD4E0CB15UL}},
    {{(int32_t)0X62F201ACUL,(int32_t)0xAECC336CUL}},
    {{(int32_t)0X41CE1E65UL,(int32_t)0x9235F2ECUL}},
    {{(int32_t)0X7641AF3DUL,(int32_t)0xCF043AB3UL}},
    {{(int32_t)0X5A82799AUL,(int32_t)0xA57D8666UL}},
    {{(int32_t)0X30FBC54DUL,(int32_t)0x89BE50C3UL}},
    {{(int32_t)0X73B5EBD1UL,(int32_t)0xC945DFECUL}},
    {{(int32_t)0X5133CC94UL,(int32_t)0x9D0DFE54UL}},
    {{(int32_t)0X1F19F97BUL,(int32_t)0x83D60412UL}},
    {{(int32_t)0X70E2CBC6UL,(int32_t)0xC3A94590UL}},
    {{(int32_t)0X471CECE7UL,(int32_t)0x9592675CUL}},
    {{(int32_t)0X0C8BD35EUL,(int32_t)0x809DC971UL}},
    {{(int32_t)0X6DCA0D14UL,(int32_t)0xBE31E19BUL}},
    {{(int32_t)0X3C56BA70UL,(int32_t)0x8F1D343AUL}},
    {{(int32_t)0XF9B82684UL,(int32_t)0x80277872UL}},
    {{(int32_t)0X6A6D98A4UL,(int32_t)0xB8E31319UL}},
    {{(int32_t)0X30FBC54DUL,(int32_t)0x89BE50C3UL}},
    {{(int32_t)0XE70747C4UL,(int32_t)0x8275A0C0UL}},
    {{(int32_t)0X66CF8120UL,(int32_t)0xB3C0200CUL}},
    {{(int32_t)0X25280C5EUL,(int32_t)0x8582FAA5UL}},
    {{(int32_t)0XD4E0CB15UL,(int32_t)0x877B7BECUL}},
    {{(int32_t)0X62F201ACUL,(int32_t)0xAECC336CUL}},
    {{(int32_t)0X18F8B83CUL,(int32_t)0x8275A0C0UL}},
    {{(int32_t)0XC3A94590UL,(int32_t)0x8F1D343AUL}},
    {{(int32_t)0X5ED77C8AUL,(int32_t)0xAA0A5B2EUL}},
    {{(int32_t)0X0C8BD35EUL,(int32_t)0x809DC971UL}},
    {{(int32_t)0XB3C0200CUL,(int32_t)0x99307EE0UL}},
    {{(int32_t)0X5A82799AUL,(int32_t)0xA57D8666UL}},
    {{(int32_t)0000000000UL,(int32_t)0x80000000UL}},
    {{(int32_t)0XA57D8666UL,(int32_t)0xA57D8666UL}},
    {{(int32_t)0X55F5A4D2UL,(int32_t)0xA1288376UL}},
    {{(int32_t)0XF3742CA2UL,(int32_t)0x809DC971UL}},
    {{(int32_t)0X99307EE0UL,(int32_t)0xB3C0200CUL}},
    {{(int32_t)0X5133CC94UL,(int32_t)0x9D0DFE54UL}},
    {{(int32_t)0XE70747C4UL,(int32_t)0x8275A0C0UL}},
    {{(int32_t)0X8F1D343AUL,(int32_t)0xC3A94590UL}},
    {{(int32_t)0X4C3FDFF4UL,(int32_t)0x99307EE0UL}},
    {{(int32_t)0XDAD7F3A2UL,(int32_t)0x8582FAA5UL}},
    {{(int32_t)0X877B7BECUL,(int32_t)0xD4E0CB15UL}},
    {{(int32_t)0X471CECE7UL,(int32_t)0x9592675CUL}},
    {{(int32_t)0XCF043AB3UL,(int32_t)0x89BE50C3UL}},
    {{(int32_t)0X8275A0C0UL,(int32_t)0xE70747C4UL}},
    {{(int32_t)0X41CE1E65UL,(int32_t)0x9235F2ECUL}},
    {{(int32_t)0XC3A94590UL,(int32_t)0x8F1D343AUL}},
    {{(int32_t)0X80277872UL,(int32_t)0xF9B82684UL}},
    {{(int32_t)0X3C56BA70UL,(int32_t)0x8F1D343AUL}},
    {{(int32_t)0XB8E31319UL,(int32_t)0x9592675CUL}},
    {{(int32_t)0X809DC971UL,(int32_t)0x0C8BD35EUL}},
    {{(int32_t)0X36BA2014UL,(int32_t)0x8C4A142FUL}},
    {{(int32_t)0XAECC336CUL,(int32_t)0x9D0DFE54UL}},
    {{(int32_t)0X83D60412UL,(int32_t)0x1F19F97BUL}},
    {{(int32_t)0X30FBC54DUL,(int32_t)0x89BE50C3UL}},
    {{(int32_t)0XA57D8666UL,(int32_t)0xA57D8666UL}},
    {{(int32_t)0X89BE50C3UL,(int32_t)0x30FBC54DUL}},
    {{(int32_t)0X2B1F34EBUL,(int32_t)0x877B7BECUL}},
    {{(int32_t)0X9D0DFE54UL,(int32_t)0xAECC336CUL}},
    {{(int32_t)0X9235F2ECUL,(int32_t)0x41CE1E65UL}},
    {{(int32_t)0X25280C5EUL,(int32_t)0x8582FAA5UL}},
    {{(int32_t)0X9592675CUL,(int32_t)0xB8E31319UL}},
    {{(int32_t)0X9D0DFE54UL,(int32_t)0x5133CC94UL}},
    {{(int32_t)0X1F19F97BUL,(int32_t)0x83D60412UL}},
    {{(int32_t)0X8F1D343AUL,(int32_t)0xC3A94590UL}},
    {{(int32_t)0XAA0A5B2EUL,(int32_t)0x5ED77C8AUL}},
    {{(int32_t)0X18F8B83CUL,(int32_t)0x8275A0C0UL}},
    {{(int32_t)0X89BE50C3UL,(int32_t)0xCF043AB3UL}},
    {{(int32_t)0XB8E31319UL,(int32_t)0x6A6D98A4UL}},
    {{(int32_t)0X12C8106FUL,(int32_t)0x8162AA04UL}},
    {{(int32_t)0X8582FAA5UL,(int32_t)0xDAD7F3A2UL}},
    {{(int32_t)0XC945DFECUL,(int32_t)0x73B5EBD1UL}},
    {{(int32_t)0X0C8BD35EUL,(int32_t)0x809DC971UL}},
    {{(int32_t)0X8275A0C0UL,(int32_t)0xE70747C4UL}},
    {{(int32_t)0XDAD7F3A2UL,(int32_t)0x7A7D055BUL}},
    {{(int32_t)0X0647D97CUL,(int32_t)0x80277872UL}},
    {{(int32_t)0X809DC971UL,(int32_t)0xF3742CA2UL}},
    {{(int32_t)0XED37EF91UL,(int32_t)0x7E9D55FCUL}},
};

static const tdct4_twd_fr32 descr = { 512, dct4_twd_512, dct3_512, rfft_256, fft_128 };
const dct_handle_t dct4_32_512=(dct_handle_t)&descr;
const dct_handle_t mdct_32_512=(dct_handle_t)&descr;
