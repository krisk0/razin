// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

/*
This file contains modified subroutine _fmpq_reconstruct_fmpz_2() owned by 
 F.Johansson 

Decreasing modulo/reusing discovered denominator trick learnt from solve1()
 subroutine implemented by V.Shoup
*/

#define LOUD_RR_IO 0
#define RR_SKIP_CHECK 1

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/longlong.h>
#include "../fmpz/fmpz_.h"
#include "fmpq_.h"

#define ROT(u,v,t)   \
    do { fmpz _t = *u; *u = *v; *v = *t; *t = _t; } while (0);

#define MulMod( _t, _d, M )  \
 {                          \
  mpz_mul( _t, _t, _d );   \
  mpz_mod( _t, _t, M );   \
 }

void __inline__ static
MulMod_fmpz(mpz_t tgt,const fmpz_t sou,const mpz_t M)
// multiply positive numbers modulo M
 {
  fmpz s=*sou;
  if( COEFF_IS_MPZ(s) )
   mpz_mul(tgt,tgt, COEFF_TO_PTR(s) );
  else
   mpz_mul_ui(tgt,tgt,s);
  mpz_mod(tgt,tgt,M);
 }

void __inline__ static
MulMod_2x(mpz_t tgt,const mpz_t sou,slong log2_M)
// multiply positive numbers modulo 2**log2_M
 {
  // TODO: implement or find a recursive algorithm to multiply gmp integers
  //  modulo degree of 2
  mpz_mul(tgt,tgt,sou);
  mpz_mod_2x(tgt,log2_M);
 }

void __inline__ static
MulMod_fmpz_2x(mpz_t tgt,const fmpz_t sou,slong log2_M)
// multiply positive numbers modulo 2**log2_M
 {
  if(fmpz_fits_si(sou))
   mpz_mul_si(tgt, tgt, fmpz_get_si(sou));
  else
   mpz_mul(tgt,tgt,COEFF_TO_PTR(*sou));
  mpz_mod_2x(tgt,log2_M);
 }

void lcm_2arg(mpz_t tgt,const fmpz_t sou)
 {
  register fmpz s=*sou;
  if(!COEFF_IS_MPZ(s))
   {
    if(s != 1)
     mpz_lcm_ui(tgt, tgt, (mp_limb_t)s);    // sou is small but not 1
   }
  else
   mpz_lcm(tgt, tgt, COEFF_TO_PTR(s));      // sou is big
 }

int __inline__ static
maybe_decrease_M(mpz_t M,mp_limb_t p,mp_limb_t log2_N,fmpz_t D)
// divide M by p so new M is still equal or greater than 2*N*D
 {
  int decreased=0;
  mpz_t b; mpz_init(b);
  fmpz_get_mpz(b,D);
  mpz_mul_2exp(b,b,1+log2_N);
  // M must be decreased iff b*p <= M
  mpz_mul_ui(b,b,p);
  // M must be decreased iff b <= M
  while( mpz_cmp(b,M) <= 0 )
   {
    decreased=1;
    mpz_divexact_ui(M,M,p);
   }
  mpz_clear(b);
  return decreased;
 }

int __inline__ static
maybe_decrease_M_2x(mpz_t Mold,slong* log2_Mold,mp_limb_t log2_N,fmpz_t D)
// possibly decrease M so it still a degree of 2 and not smaller than 2*N*D
 {
  slong log2_M=fmpz_clog_ui(D,2)+1+log2_N;
  if(log2_M < *log2_Mold)
   {
    mpz_set_ui(Mold,1); mpz_mul_2exp(Mold,Mold,log2_M);
    *log2_Mold = log2_M;
    return 1;
   }
  return 0;
 }

int __inline__ static
fmpz_cmpabs_log2(const fmpz_t x,mp_limb_t y)
 {
  register slong xx=(slong)(*x);
  if(!COEFF_IS_MPZ(xx))
   return -1;          // |x|<2**64, 2**y>2**64
  mpz_ptr z=COEFF_TO_PTR(xx);
  xx=mpz_size(z);      // xx>0
  if(xx*FLINT_BITS <= y)
   return -1;          // x is too short, no need to further compare
  --xx;
  if( xx*FLINT_BITS > y )
   return 1;           // even if senior limb is 1, x is larger
  int lz;
  count_leading_zeros(lz,mpz_getlimbn(z,xx));
  xx=(xx*FLINT_BITS)+FLINT_BITS-lz;   // xx=|x| bit_length
  if(xx <= y)
   return -1;     // |x| fits y bits
  if(xx-1 > y)
   return 1;      // |x| > 2**y because its bit-length is y+2 or more
  // Now either |x| is degree of 2, and result should be 0
  //  or |x| is not degree of 2, and result should be 1
  return !abs_x_is_degree_of_2(z);
 }

//fmpz_cmp_log2 is more expensive than fmpz_cmpabs_log2
int __inline__ static
fmpz_cmp_log2(const fmpz_t x,mp_limb_t y)
 {
  if(fmpz_cmp_ui(x,0)<=0)
   return -1;                    // x<=0, 2**y>0
  return fmpz_cmpabs_log2(x,y);
 }

// replace fmpz_cmp with fmpz_cmpabs because a is non-negative
#define RR_tome0                   \
    fmpz_t q, r, s, t;              \
    int success = 0;                  \
    if (fmpz_cmpabs_log2(a, N2) <= 0)   \
    {                                 \
        fmpz_set(n, a);              \
        fmpz_one(d);                 \
        return 1;                    \
    }                                \
    fmpz_sub(n, a, m);                \
    if (fmpz_cmpabs_log2(n, N2) <= 0) \
    {                                 \
        fmpz_one(d);                 \
        return 1;                    \
    }                                \
    fmpz_init(q);                     \
    fmpz_init(r);                       \
    fmpz_init(s);                          \
    fmpz_init(t);                             \
    fmpz_set(r, m); fmpz_zero(s);                \
    fmpz_set(n, a); fmpz_one(d);                    \
    while (fmpz_cmpabs_log2(n, N2) > 0)               \
    {                                                   \
        fmpz_fdiv_q(q, r, n);                             \
        fmpz_mul(t, q, n); fmpz_sub(t, r, t); ROT(r, n, t); \
        fmpz_mul(t, q, d); fmpz_sub(t, s, t); ROT(s, d, t); \
    }                                                      \
    if (fmpz_sgn(d) < 0)                                  \
    {                                                    \
        fmpz_neg(n, n);                                 \
        fmpz_neg(d, d);                                \
    }

#define RR_tome1     \
    fmpz_clear(q);   \
    fmpz_clear(r);   \
    fmpz_clear(s);   \
    fmpz_clear(t);   \
    return success;

int reconstruct_rational_log2_plain(fmpz_t n, fmpz_t d,
  const fmpz_t a, const fmpz_t m, mp_limb_t N2, const fmpz_t D,
  int skip_check)
/*
like _fmpq_reconstruct_fmpz_2(), but 
 * numerator bound is logarithmic
 * optionally skip check
 * dont run gcd algorithm when performing check 
*/
 {
  RR_tome0
  if (fmpz_cmp(d, D) <= 0)
   {
    if(skip_check)
     success=1;
    else
     {
      /* 
       a=n/d modulo m  iff gcd(n,d)=1. But it should be faster to directly
       compare a*d to n
      */
      fmpz_mul(s, a,d);
      fmpz_mod(t, s,m);
      if(fmpz_cmp_ui(n,0)>0)
       success = !fmpz_cmp(n,t);
      else
       {
        fmpz_mod(s, n,m);
        success = !fmpz_cmp(s,t);
       }
     }
   }
  RR_tome1
 }

int reconstruct_rational_log2_log2(fmpz_t n, fmpz_t d,
  const fmpz_t a, const fmpz_t m, mp_limb_t N2, mp_limb_t D2,
  int skip_check)
// like reconstruct_rational_log2_plain(), but denominator bound is logarithmic
 {
  RR_tome0
  #if LOUD_RR_IO
   fmpz_hex_print("RR inside n=",n,1);
   fmpz_hex_print("RR inside d=",d,1);
   flint_printf("d in range: %d\n",fmpz_cmpabs_log2(d, D2));
  #endif
  if (fmpz_cmpabs_log2(d, D2) <= 0)  // d is positive, using cheaper
   {                                 //  fmpz_cmpabs_log2
    if(skip_check)
     success=1;
    else
     {
      fmpz_mul(s, a,d);
      fmpz_mod(t, s,m);
      #if LOUD_RR_IO
       fmpz_hex_print("          a*d=",t,1);
      #endif
      if(fmpz_cmp_ui(n,0)>0)
       success = !fmpz_cmp(n,t);
      else
       {
        fmpz_mod(s, n,m);
        #if LOUD_RR_IO
         fmpz_hex_print("          n mod M=",s,1);
        #endif
        success = !fmpz_cmp(s,t);
       }
     }
   }
  RR_tome1
 }

#undef RR_tome1
#undef RR_tome0

int __inline__ static
reconstruct_rational_take_denominator(fmpz_t d,mpz_t aa,const mpz_t mm,
  mp_limb_t log2_N, fmpz_t D, int skip_check)
// perform rational reconstruction, drop numerator, return denominator
 {
  fmpz_t m; fmpz_init_set_readonly(m, mm); 
  fmpz_t a; fmpz_init_set_readonly(a, aa); 
  fmpz_t n; fmpz_init(n);
  int r=reconstruct_rational_log2_plain(n,d, a,m, log2_N,D, skip_check);
  #if LOUD_RR_IO
   flint_printf("RR log2(N)=%llX\n",log2_N);
   fmpz_hex_print("RR D=",D,1);
   fmpz_hex_print("RR a=",a,1);
   fmpz_hex_print("RR M=",m,1);
   if(r)
    {
     fmpz_hex_print("   n=",n,1);
     fmpz_hex_print("   d=",d,1);
    }
   else
    flint_printf("RR failed\n");
   flint_printf("\n");
  #endif
  fmpz_clear(n);
  fmpz_clear_readonly(a);
  fmpz_clear_readonly(m);
  return r;
 }

int __inline__ static
reconstruct_rational_take_denominator_log2(fmpz_t d,mpz_t aa,const mpz_t mm,
  mp_limb_t log2_N, mp_limb_t log2_D, int skip_check)
// perform rational reconstruction, drop numerator, return denominator
 {
  fmpz_t m; fmpz_init_set_readonly(m, mm); 
  fmpz_t a; fmpz_init_set_readonly(a, aa); 
  fmpz_t n; fmpz_init(n);
  int r=reconstruct_rational_log2_log2(n,d, a,m, log2_N,log2_D, skip_check);
  #if LOUD_RR_IO
   flint_printf("RR log2(N)=%llX\n",log2_N);
   flint_printf("RR log2(D)=%llX\n",log2_D);
   fmpz_hex_print("RR a=",a,1);
   fmpz_hex_print("RR M=",m,1);
   if(r)
    {
     fmpz_hex_print("   n=",n,1);
     fmpz_hex_print("   d=",d,1);
    }
   else
    flint_printf("RR failed\n");
   flint_printf("\n");
  #endif
  fmpz_clear(n);
  fmpz_clear_readonly(a);
  fmpz_clear_readonly(m);
  return r;
 }

void
det_divisor_rational_reconstruction(mpz_t d,mpz_ptr x,mpz_t M,mp_limb_t p,
  slong n,mp_limb_t log2_N, mp_limb_t log2_D)
/*
M=modulo=p**k, for some integer k
log2_N: log2(upper bound on numerators)
log2_D: log2(upper bound on denominators)
x[i] is in range 0..M-1, for i in 0..n-1

Algorithm behind this subroutine inspired by Victor Shoup solve1() subroutine 
 found in mat_ZZ.c (which is part of NTL)
*/
 {
  #if RAT_REC_TAKES_D_SERIOUSLY==0
   mpz_set_ui(d,1);
  #endif
  mpz_t d_mod_M; mpz_init_set_ui(d_mod_M,1);
  slong i;
  int M_modified=0,rc;
  mpz_ptr xI;
  fmpz_t found; fmpz_init(found);
  fmpz_t D; fmpz_init(D); // D=0
  #if RAT_REC_TAKES_D_SERIOUSLY
   // if d>1 then mix it into x[]
   if( mpz_cmp_ui(d,1)>0 )
    {
     MulMod( x, d, M );
    }
  #endif
  for(i=0,xI=x; i<n; i++,xI++)
   {
    if(M_modified)
     mpz_mod(xI,xI,M);
    if( mpz_cmp_ui(d_mod_M,1) )
     MulMod( xI, d_mod_M, M );
    if( fmpz_cmp_ui(D,0) )
     rc=reconstruct_rational_take_denominator(found, xI, M, log2_N, D, 
      RR_SKIP_CHECK);
    else
     rc=reconstruct_rational_take_denominator_log2(found, xI, M, log2_N, 
      log2_D, RR_SKIP_CHECK);
    if(!rc)
     {
      flint_printf("Exception (det_divisor_rational_reconstruction): "
                     "Rational reconstruction failed.\n");
      abort();
     }
    if( fmpz_cmp_ui(found,1) )
     {
      mul_mpz_fmpz( d, found );
      if( i != n-1 )
       {
        // Decrease denominator bound, maybe decrease M
        if( !fmpz_cmp_ui(D,0) )
         {
          fmpz_set_ui(D,1); fmpz_mul_2exp(D,D,log2_D);
         }
        fmpz_fdiv_q(D, D, found);
        M_modified |= maybe_decrease_M(M,p,log2_N,D);
        // multiply d_mod_M by found and reduce result modulo M
        if(M_modified)
         {
          mpz_mod(d_mod_M,d_mod_M,M);
          fmpz_mod_2arg(found,M);
         }
        MulMod_fmpz(d_mod_M, found, M);
       }
     }
   }
  fmpz_clear(D);
  fmpz_clear(found);
  mpz_clear(d_mod_M);
 }

void
rational_reconstruction_2deg(mpz_t d,mpz_ptr x,slong n,mpz_t M,slong log2_M,
  mp_limb_t log2_N,mp_limb_t log2_D)
/*
d: on entry = 1, on exit common denominator of reconstructed salvation
x: vector of length n, x*A equals B modulo M, 0 <= x[i] < M, det A is odd
M=2**log2_M
log2_N: upper approximation to log2(numerator bound)
log2_D: upper approximation to log2(denominator bound)
 
essentially the same algorithm as det_divisor_rational_reconstruction() 
*/
 {
  mpz_t d_mod_M; mpz_init_set_ui(d_mod_M,1);
  fmpz_t found; fmpz_init(found);
  fmpz_t D; fmpz_init(D); // D=0
  slong i;
  int M_modified=0,rc;
  mpz_ptr xI;
  for(i=0,xI=x; i<n; i++,xI++)
   {
    if(M_modified)
     mpz_mod_2x(xI,log2_M);
    if( mpz_cmp_ui(d_mod_M,1) )
     MulMod_2x( xI, d_mod_M, log2_M );
    if( fmpz_cmp_ui(D,0) )
     rc=reconstruct_rational_take_denominator(found, xI, M, log2_N, D, 
      RR_SKIP_CHECK);
    else
     rc=reconstruct_rational_take_denominator_log2(found, xI, M, log2_N, 
      log2_D, RR_SKIP_CHECK);
    if(!rc)
     {
      flint_printf("Exception (det_divisor_rational_reconstruction): "
                     "Rational reconstruction failed.\n");
      abort();
     }
    if( fmpz_cmp_ui(found,1) )
     {
      mul_mpz_fmpz( d, found );
      if( i != n-1 )
       {
        // Decrease denominator bound, maybe decrease M
        if( !fmpz_cmp_ui(D,0) )
         {
          fmpz_set_ui(D,1); fmpz_mul_2exp(D,D,log2_D);
         }
        fmpz_fdiv_q(D, D, found);
        M_modified |= maybe_decrease_M_2x(M,&log2_M,log2_N,D);
        // multiply d_mod_M by found and reduce result modulo M
        if(M_modified)
         {
          mpz_mod_2x(d_mod_M,log2_M);
          fmpz_mod_2x(found,log2_M);
         }
        MulMod_fmpz_2x(d_mod_M, found, log2_M);
       }
     }
   }
  fmpz_clear(D);
  fmpz_clear(found);
  mpz_clear(d_mod_M);
 }

#undef ROT
#undef RR_SKIP_CHECK
#undef MulMod
#undef NDEBUG
#undef LOUD_RR_IO
