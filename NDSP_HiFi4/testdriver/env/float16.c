/* ------------------------------------------------------------------------ */
/* Copyright (c) 2018 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2018 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
 * Half precision floating-point arithmetic and manipulation functions.
 */

/* Partable data types. */
#include "types.h"
/* Half precision floating-point arithmetic and manipulation functions. */
#include "float16.h"

#define MIN(a,b)  ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)  ( (a)>(b) ? (a) : (b) )

static int16_t _S_exp0_l (int32_t x)
{
  int16_t z=0;
  if ( x==0 )  return 0x1F;
  while ( (int32_t)(x^(x<<1))>0 ) //MSB != MSB-1
  {
    x<<=1;
    z++;
  }
  return z;
}

/* Convert single precision floating-point data to half precision format. */
float16_t conv_f32_to_f16( float32_t x )
{
  union ufloat32uint32 _x;
  int32_t s, e, f;
  int32_t tb, rbit, tail;

  _x.f = x;

  s = (int32_t)(_x.u>>31);
  e = ( (int32_t)(_x.u>>23) & 255L ) - 127L;
  
  /* Restore the implicit leading significand bit. */
  f = (int32_t)( _x.u & ( (1L<<23) - 1 ) );

  if ( 128==e && 0!=f ) /* NaN */
  {
    /* Return a quiet (!) NaN with original payload bits. */
    return (float16_t)( (s<<15) | 0x7e00 | (f>>(23-10)) ) ; 
  }
  else if ( e>15 ) /* Too large magnitude */
  {
    return (float16_t)( (s<<15) | 0x7c00 );
  }
  else if ( e<-25 ) /* Too small magnitude */
  {
    return (float16_t)(s<<15);
  }
  else if ( -14>e ) /* Subnormal half precision value */
  {
    /* Restore the implicit leading significand bit. */ 
    f += (1L<<23);
    /* 10 significand trailing bits */
    tb = ( f >> (23-10-e-14) ); 
    /* Rounding bit value */
    rbit = ( f >> (23-10-1-e-14) ) & 1L; 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1L<<(23-10-1-e-14)) - 1 ) ) );
    e = -15;
  }
  else /* Normal half precision value */
  {
    /* 10 significand trailing bits */
    tb = ( f >> (23-10) ); 
    /* Rounding bit value */
    rbit = ( f >> (23-10-1) ) & 1L; 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1L<<(23-10-1)) - 1 ) ) );
  }

  if ( 0!=rbit ) 
  {
    /* Round to nearest. In case of a tie round to even mantissa. */
    tb += ( 0==tail ? (tb&1) : 1 );
    /* Look if the exponent shall be adjusted for the rounded mantissa. */
    if ( tb & (1L<<10) )
    { 
      /* Correct the exponent and mantissa, and check if the rounded value overflows. */
      tb = ( ++e<=15 ? ( (tb+1024) >> 1 ) : 0 );
    }
  }

  return (float16_t)( (s<<15) | ( (e+15) << 10 ) | ( tb & 0x3ff ) );

} /* conv_f32_to_f16() */

/* Convert double precision floating-point data to half precision format. */
float16_t conv_f64_to_f16( float64_t x )
{
  union ufloat64uint64 _x;
  int32_t s, e;
  uint64_t f;
  int32_t tb, rbit, tail;

  _x.f = x;

  s = (int32_t)(_x.u>>63);
  e = ( (int32_t)(_x.u>>52) & 2047L ) - 1023L;
  f = ( _x.u & ( (1ULL<<52) - 1 ) );

  if ( 1024==e && 0!=f ) /* NaN */
  {
    /* Return a quiet (!) NaN with original payload bits. */
    return (float16_t)( (s<<15) | 0x7e00 | (int32_t)(f>>(52-10)) ) ; 
  }
  else if ( e>15 ) /* Too large magnitude */
  {
    return (float16_t)( (s<<15) | 0x7c00 );
  }
  else if ( e<-25 ) /* Too small magnitude */
  {
    return (float16_t)(s<<15);
  }
  else if ( -14>e ) /* Subnormal half precision value */
  {
    /* Restore the implicit leading significand bit. */ 
    f += (1ULL<<52);
    /* 10 significand trailing bits */
    tb = (int32_t)(f>>(52-10-e-14)); 
    /* Rounding bit value */
    rbit = (int32_t)( (f>>(52-10-1-e-14)) & 1L ); 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1ULL<<(52-10-1-e-14)) - 1 ) ) );
    e = -15;
  }
  else /* Normal half precision value */
  {
    /* 10 significand trailing bits */
    tb = (int32_t)(f>>(52-10)); 
    /* Rounding bit value */
    rbit = (int32_t)( (f>>(52-10-1)) & 1L ); 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1ULL<<(52-10-1)) - 1 ) ) );
  }

  if ( 0!=rbit ) 
  {
    /* Round to nearest. In case of a tie round to even mantissa. */
    tb += ( 0==tail ? (tb&1) : 1 );
    /* Look if the exponent shall be adjusted for the rounded mantissa. */
    if ( tb & (1L<<10) )
    { 
      /* Correct the exponent and mantissa, and check if the rounded value overflows. */
      tb = ( ++e<=15 ? ( (tb+1024) >> 1 ) : 0 );
    }
  }

  return (float16_t)( (s<<15) | ( (e+15) << 10 ) | ( tb & 0x3ff ) );

} /* conv_f64_to_f16() */

/* Convert half precision floating-point data to single precision format. */
float32_t conv_f16_to_f32( float16_t x )
{
  int32_t s,e,f;
  union ufloat32uint32 y = {0};

  s = ( (uint16_t)x >> 15 );
  e = ( (uint16_t)x >> 10 ) & 31;
  f = ( (uint16_t)x >>  0 ) & 1023;

  if ( 31==e )
  {
    if ( 0==f ) /* +/-Infinity */
    {
      y.u = ( (uint32_t)s << 31 ) | ( 0xffUL << 23 );
    }
    else /* NaN */
    {
      /* Return a quiet (!) NaN with original payload bits. */
      y.u = ( (uint32_t)s << 31 ) | ( 0x1ffUL << 22 ) | ( ( (uint32_t)f & 0x1ffUL ) << (23-10) ); 
    }
  }
  else if ( 0==e && 0!=f ) /* Subnormal */
  {
    int nsa = _S_exp0_l( f );
    /* Clear the leading bit of the significand. */
    f = (f<<(nsa-7)) ^ (1L<<23);
    /* Q(23-e) <- [ Q(10+14) + nsa ] + (23-30) */
    y.u = ( (uint32_t)s << 31 ) | ( (uint32_t)(6-nsa+127) << 23 ) | (uint32_t)f;
  }
  else if ( 0==e ) /* +/-0 */
  {
    y.u = ( (uint32_t)s << 31 );
  }
  else /* Normal finite value. */
  {
    y.u = ( (uint32_t)s << 31 ) | ( (uint32_t)(e-15+127) << 23 ) | ( (uint32_t)f << (23-10) );
  }

  return (y.f);

} /* conv_f16_to_f32() */

/* Convert half precision floating-point data to double precision format. */
float64_t conv_f16_to_f64( float16_t x )
{
  int32_t s,e,f;
  union ufloat64uint64 y = {0};

  s = ( (uint16_t)x >> 15 );
  e = ( (uint16_t)x >> 10 ) & 31;
  f = ( (uint16_t)x >>  0 ) & 1023;

  if ( 31==e )
  {
    if ( 0==f ) /* +/-Infinity */
    {
      y.u = ( (uint64_t)s << 63 ) | ( 0x7ffULL << 52 );
    }
    else /* NaN */
    {
      /* Return a quiet (!) NaN with original payload bits. */
      y.u = ( (uint64_t)s << 63 ) | ( 0xfffULL << 51 ) | ( ( (uint64_t)f & 0x1ffULL ) << (52-10) ); 
    }
  }
  else if ( 0==e && 0!=f ) /* Subnormal */
  {
    int nsa = _S_exp0_l( f );
    /* Clear the leading bit of the significand. */
    f = (f<<nsa) & ~(1L<<30);
    /* Q(52-e) <- [ Q(10+14) + nsa ] + (52-30) */
    y.u = ( (uint64_t)s << 63 ) | ( (uint64_t)(6-nsa+1023) << 52 ) | ( (uint64_t)f << (52-30) );
  }
  else if ( 0==e ) /* +/-0 */
  {
    y.u = ( (uint64_t)s << 63 );
  }
  else /* Normal finite value. */
  {
    y.u = ( (uint64_t)s << 63 ) | ( (uint64_t)(e-15+1023) << 52 ) | ( (uint64_t)f << (52-10) );
  }

  return (y.f);

} /* conv_f16_to_f64() */

/* Sum of two half precision floating-point values */
float16_t add_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) +
                            conv_f16_to_f64(y) ) );
}

/* Sum of two half precision floating-point values */
float16_t sub_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) -
                            conv_f16_to_f64(y) ) );
}

/* Multiply two half precision floating-point values */
float16_t mul_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) *
                            conv_f16_to_f64(y) ) );
}

/* Multiply two half precision floating-point values */
float16_t div_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) /
                            conv_f16_to_f64(y) ) );
}

/* Fused multiply-add: return x*y+z. */
float16_t fma_f16( float16_t x, float16_t y, float16_t z )
{
  float64_t xd, yd, zd;

  xd = conv_f16_to_f64(x);
  yd = conv_f16_to_f64(y);
  zd = conv_f16_to_f64(z);

  return ( conv_f64_to_f16( xd*yd+zd ) );
}

/* Return absolute value of half precision floating-point number. */
float16_t fabs_f16( float16_t x )
{
  return ( isnan_f16(x) ? x : (float16_t)( (uint16_t)x & 0x7fff ) );
}

/* Multiply the first argument by two, raised to the power of the second argument. */
float16_t ldexp_f16( float16_t x, int n )
{
  int n0, n1, n2;
  float16_t y, s0, s1, s2;

  /* Sensible range of exponent adjustment is [-41,40], as shown below:
   *  - max finite -> zero: pow2(it_half((2-2^-10)*2^15),-41) == 0
   *  - min subnormal -> Infinity: pow2(it_half(2^-24),40) == Infinity */
  n = MAX( -41, MIN( 40, n ) );

  /* Multiplication by a power of 2 is the most convenient way to modify
   * the exponent. Maximum representable power of two for the half precision
   * is 15, thus we have to split n into 3 components to cover the full range
   * of exponent adjustment: [-41,40]. The value is partitioned in such a way
   * that a subnormal result may appear only once in three products, excepting
   * trivial multiplucations by 1. This guarantees that the subnormal result is
   * rounded only once. */
  n2 = MAX( -14, MIN( 15, n    ) ); s2 = (float16_t)( (uint16_t)(n2+15)<<10 );
  n1 = MAX( -14, MIN( 15, n-n2 ) ); s1 = (float16_t)( (uint16_t)(n1+15)<<10 );
  n0 = n-n2-n1;                     s0 = (float16_t)( (uint16_t)(n0+15)<<10 );

  /* 1. The order of multiplications is crucial.
   * 2. We rely on the * operator to raise the FE_OVERFLOW floating-point
   *    exception whenever result overflows. */
  y = mul_f16( mul_f16( mul_f16(x,s0), s1 ), s2 );

  return (y);

} /* ldexp_f16() */

typedef enum { _lt_f16, _lte_f16, _gt_f16, 
               _gte_f16, _eq_f16, _oneq_f16 } f16_order_t;

static int compare_f16( float16_t a, float16_t b, f16_order_t ord )
{
  int r;
  int32_t ai, bi;
  if ( un_f16(a,b) ) return (0);
  ai = ( (uint16_t)a & 0x7fff ); if ( (uint16_t)a & 0x8000 ) ai = -ai;
  bi = ( (uint16_t)b & 0x7fff ); if ( (uint16_t)b & 0x8000 ) bi = -bi;
  r = ( ( ord == _lt_f16 ) ? (ai<bi ) : ( ord == _lte_f16  ) ? (ai<=bi) :
        ( ord == _gt_f16 ) ? (ai>bi ) : ( ord == _gte_f16  ) ? (ai>=bi) :
        ( ord == _eq_f16 ) ? (ai==bi) : ( ord == _oneq_f16 ) ? (ai!=bi) : -1 );
  NASSERT( r>=0 );
  return (r);
}

int lt_f16  ( float16_t a, float16_t b ) { return compare_f16( a, b, _lt_f16   ); } /* a<b  */
int lte_f16 ( float16_t a, float16_t b ) { return compare_f16( a, b, _lte_f16  ); } /* a<b  */
int gt_f16  ( float16_t a, float16_t b ) { return compare_f16( a, b, _gt_f16   ); } /* a<b  */
int gte_f16 ( float16_t a, float16_t b ) { return compare_f16( a, b, _gte_f16  ); } /* a<b  */
int eq_f16  ( float16_t a, float16_t b ) { return compare_f16( a, b, _eq_f16   ); } /* a<b  */
int oneq_f16( float16_t a, float16_t b ) { return compare_f16( a, b, _oneq_f16 ); } /* a!=b */
