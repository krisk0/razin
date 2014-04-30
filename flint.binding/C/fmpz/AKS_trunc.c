// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// 
// Licence: public domain

// This file is based on public domain code put online by 袁轶君 23 apr 2014 at 
//  https://github.com/YijunYuan/AKS_FLINT

// Original code did not work for large n, set `c` border on stage 2 and 5 
//  incorrectly. So I (Денис Крыськов) rewrote it

#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fmpz_mod_poly.h>
#include<flint/arith.h>
//#include<mpir.h>     gmp or mpir header is included by FLINT
//#include<stdio.h>    stdio.h too
#include<math.h>

int 
multiplicative_order_greater( mp_limb_t r, mp_limb_t n, mp_limb_t c)
/*
return 1 iff n belongs to (Z/rZ)* and its order is
 equal or greater than c

c > 1

r < ceil( log2(n)**5 ) <= 64**5 = 2*30, so it is safe to multiply then reduce
*/
 {
  mp_limb_t m=n % r;
  if( n_gcd(r,m) != 1 )
   return 0;
  n=m;
  mp_limb_t deg=1;
  while(1)
   {
    // n**deg equals m modulo r, deg < c
    if(m==1)
     return 0;
    if( ++deg >= c )
     return 1;
    m = (m*n) % r;
   }
 }

mp_limb_t
ceil_square_log2(mp_limb_t n)
// count ceil( log2(n) ** 2 )
 {
  slong prec=50;
  mpfr_t lo,up; mpfr_init( lo ); mpfr_init( up );
  mp_limb_t r;
  while(1)
   {
    mpfr_set_ui( lo, n,    MPFR_RNDD);  mpfr_set_ui( up, n,    MPFR_RNDU);
    mpfr_log2( lo, lo,     MPFR_RNDD);  mpfr_log2( up, up,     MPFR_RNDU);
    mpfr_mul ( lo, lo, lo, MPFR_RNDD);  mpfr_mul ( up, up, up, MPFR_RNDU);
    mpfr_ceil( lo, lo );                mpfr_ceil( up, up );
    r=mpfr_get_ui(lo, MPFR_RNDN);
    if( r == mpfr_get_ui(up, MPFR_RNDN) )
     {
      mpfr_clear(lo); mpfr_clear(up);
      return r;
     }
    prec += 50;
    if( prec > MPFR_PREC_MAX )
     {
      printf("problem in ceil_square_log2()\n"); 
      abort();
     }
    mpfr_set_prec( lo, prec ); mpfr_set_prec( up, prec );
   }
 }

mp_limb_t
stage5_AKS_bound(mp_limb_t r, mp_limb_t n)
// count floor( sqrt(phi(r)) * log2(n) )
 {
  mp_limb_t k=n_euler_phi( r );
  mpfr_t lo0,up0; mpfr_init( lo0 ); mpfr_init( up0 );
  mpfr_t lo1,up1; mpfr_init( lo1 ); mpfr_init( up1 );
  slong prec=50;
  mp_limb_t q;
  while(1)
   {
    mpfr_sqrt_ui( lo0, k,    MPFR_RNDD);  mpfr_sqrt_ui( up0,   k,   MPFR_RNDU);
    mpfr_set_ui(  lo1, n,    MPFR_RNDD);  mpfr_set_ui(  up1,   n,   MPFR_RNDU);
    mpfr_log2(    lo1, lo1,  MPFR_RNDD);  mpfr_log2(    up1,  up1,  MPFR_RNDU);
    mpfr_mul( lo0, lo0, lo1, MPFR_RNDD);  mpfr_mul( up0, up0, up1, MPFR_RNDU);
    mpfr_floor( lo0, lo0 );               mpfr_floor( up0, up0 );
    q=mpfr_get_ui(lo0, MPFR_RNDN);
    if( q == mpfr_get_ui(up0, MPFR_RNDN) )
     {
      mpfr_clear(lo0); mpfr_clear(up0);
      mpfr_clear(lo1); mpfr_clear(up1);
      return q;
     }
    prec += 50;
    if( prec > MPFR_PREC_MAX )
     {
      printf("problem in stage5_AKS_bound()\n"); 
      abort();
     }
    mpfr_set_prec( lo0, prec ); mpfr_set_prec( up0, prec );
    mpfr_set_prec( lo1, prec ); mpfr_set_prec( up1, prec );
   }
 }

int
AKS(fmpz_t n)
/*
implementation of algorithm described in "PRIMES is in P" by Manindra Agrawal
 Neeraj Kayal, Nitin Saxena

TOO SLOW for big numbers; might run out of memory due to big degree polynoms
 
n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64
*/ 
{
 /*Step1*/
 mpz_t k;mpz_init(k);fmpz_get_mpz(k,n);
 if(mpz_perfect_power_p(k)!=0)
  {
   mpz_clear(k);
   return 0;
  }
 mpz_clear(k);
 /*Step2*/
 mp_limb_t n_ui=fmpz_get_ui(n);
 mp_limb_t c_ui=ceil_square_log2(n_ui); // c_ui=(log2(n))**2 rounded up, c in 2..2**12
 mp_limb_t r_ui;
 for(r_ui=c_ui+1;;r_ui++)
  {
   if( multiplicative_order_greater( r_ui, n_ui, c_ui ) )
    break;
  }
 /*Step3*/
 // r_ui < 2**30
 mp_limb_t a_ui;
 for( a_ui=3; a_ui<=r_ui; a_ui++)
  {
   c_ui=n_gcd_full( a_ui, n_ui );
   if( (c_ui > 1) && (c_ui < n_ui) )
    return 0;
  }
 /*Step4*/
 if( n_ui <= r_ui )
  return 1;
 /*Step5*/
 c_ui = stage5_AKS_bound( r_ui, n_ui );
 fmpz_mod_poly_t modulo; fmpz_mod_poly_init(modulo,n);
 fmpz_mod_poly_set_coeff_ui(modulo,0,n_ui-1);
 fmpz_mod_poly_set_coeff_ui(modulo,r_ui, 1);  //modulo=x^r-1
 fmpz_mod_poly_t p     ; fmpz_mod_poly_init(p     ,n);
 fmpz_mod_poly_t q     ; fmpz_mod_poly_init(q     ,n);
 mp_limb_t n_modulo_r=n_ui % r_ui;
 // TODO: move more code outside loop. I won't do that since this the code is
 //  too slow
 for( a_ui=c_ui; a_ui--; ) 
  {
   fmpz_mod_poly_zero(p);
   fmpz_mod_poly_set_coeff_ui(p, 0 , a_ui);
   fmpz_mod_poly_set_coeff_ui(p, 1 , 1);     //p=x+a
   fmpz_mod_poly_zero(q);
   fmpz_mod_poly_set_coeff_ui( q, n_modulo_r, 1     );
   fmpz_mod_poly_set_coeff_ui( q,          0, a_ui  );
   fmpz_mod_poly_powmod_ui_binexp(p,p,n_ui,modulo);//p=p^n mod modulo
   if(fmpz_mod_poly_equal(p,q)==0)
    {
     fmpz_mod_poly_clear(modulo); 
     fmpz_mod_poly_clear(p); 
     fmpz_mod_poly_clear(q);
     return 0;
    }
  }
 fmpz_mod_poly_clear(modulo); 
 fmpz_mod_poly_clear(p); 
 fmpz_mod_poly_clear(q);
 return 1;
}

int 
AKS_ui(mp_limb_t n)
// n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64
 {
  fmpz_t m; fmpz_init_set_ui( m, n );
  int rc=AKS(m);
  fmpz_clear(m);
  return rc;
 }
