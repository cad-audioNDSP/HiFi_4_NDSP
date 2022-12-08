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
#include "scl_sqrt_table.h"
#include "common.h"
#ifdef COMPILER_MSVC

/*
MATLAB code:
ord=128;
xt=0:(1/ord):1;
for (i=1:(ord))
    xtbl(i)=xt(i)+((xt(i+1)-xt(i))/2);
end;
ytbl=round(sqrt(xtbl)*2^31);
ctbl=round(1./(2*sqrt(xtbl))*2^31);
ytbl=max(-2^31,min(2^31-1,ytbl));
ctbl=max(-2^31,min(2^31-1,ctbl));
*/
const int32_t _declspec(align(8)) sqrt_table[] =
#else
const int32_t                     sqrt_table[] __attribute__((aligned(8))) =
#endif
{
 134217728,	2147483647,	 232471924, 	2147483647,		
 300119964,	2147483647,	 355106730,	  2147483647,		
 402653184,	2147483647,	 445149844, 	2147483647,		
 483928900,	2147483647,	 519823025,	  2147483647,		
 553393869,	2147483647,	 585041513, 	2147483647,		
 615062898,	2147483647,	 643685611,	  2147483647,		
 671088640,	2147483647,	 697415773, 	2147483647,		
 722784585,	2147483647,	 747292683, 	2147483647,		
 771022147,	2147483647,	 794042787, 	2147483647,		
 816414567,	2147483647,	 838189443, 	2147483647,		
 859412787,	2147483647,	 880124500, 	2147483647,		
 900359891,	2147483647,	 920150384, 	2147483647,		
 939524096,	2147483647,	 958506298, 	2147483647,		
 977119809,	2147483647,	 995385311, 	2147483647,		
1013321625,	2147483647,	1030945931,	  2147483647,		
1048273967,	2147483647,	1065320189,	  2147483647,		
1082097918,	2130900515,	1098619452,	  2098855072,		
1114896182,	2068213208,	1130938678,	  2038875364,		
1146756771,	2010751598,	1162359621,	  1983760420,		
1177755783,	1957827796,	1192953261,	  1932886296,		
1207959552,	1908874354,	1222781696,	  1885735628,		
1237426310,	1863418444,	1251899625,	  1841875310,		
1266207514,	1821062491,	1280355523,	  1800939636,		
1294348895,	1781469447,	1308192592,	  1762617387,		
1321891318,	1744351429,	1335449532,	  1726641819,		
1348871473,	1709460876,	1362161168,	  1692782810,		
1375322451,	1676583559,	1388358974,	  1660840642,		
1401274219,	1645533028,	1414071510,	  1630641020,		
1426754019,	1616146146,	1439324782,	  1602031062,		
1451786701,	1588279468,	1464142555,	  1574876026,		
1476395008,	1561806289,	1488546612,	  1549056637,		
1500599818,	1536614214,	1512556978,	  1524466875,		
1524420351,	1512603139,	1536192112,	  1501012140,		
1547874349,	1489683584,	1559469076,	  1478607716,		
1570978229,	1467775280,	1582403676,	  1457177486,		
1593747216,	1446805984,	1605010588,	  1436652834,		
1616195466,	1426710480,	1627303469,	  1416971728,		
1638336161,	1407429723,	1649295054,	  1398077927,		
1660181608,	1388910104,	1670997238,	  1379920300,		
1681743312,	1371102827,	1692421154,	  1362452250,		
1703032049,	1353963368,	1713577240,	  1345631207,		
1724057932,	1337451002,	1734475296,	  1329418191,		
1744830464,	1321528399,	1755124538,	  1313777432,		
1765358587,	1306161267,	1775533649,	  1298676040,		
1785650732,	1291318043,	1795710816,	  1284083712,		
1805714853,	1276969620,	1815663770,	  1269972473,		
1825558469,	1263089103,	1835399826,	  1256316458,		
1845188694,	1249651603,	1854925906,	  1243091706,		
1864612269,	1236634043,	1874248572,	  1230275986,		
1883835584,	1224014999,	1893374053,	  1217848637,		
1902864709,	1211774541,	1912308264,	  1205790433,		
1921705413,	1199894112,	1931056833,	  1194083452,		
1940363185,	1188356400,	1949625114,	  1182710970,		
1958843251,	1177145240,	1968018211,	  1171657354,		
1977150595,	1166245512,	1986240991,	  1160907976,		
1995289972,	1155643060,	2004298098,	  1150449133,		
2013265920,	1145324612,	2022193972,	  1140267967,		
2031082780,	1135277711,	2039932856,	  1130352405,		
2048744702,	1125490652,	2057518809,	  1120691096,		
2066255659,	1115952423,	2074955721,	  1111273357,		
2083619457,	1106652658,	2092247318,	  1102089122,		
2100839745,	1097581581,	2109397173,	  1093128899,		
2117920024,	1088729972,	2126408716,	  1084383727,		
2134863654,	1080089122,	2143285240,	  1075845140		
};
