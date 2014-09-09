// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// This program tests function cmp_positive_log2() against macro fmpz_cmpabs()

// To compile: gcc -O2 -march=native -lflint -lgmp

#include "C/fmpz/fmpz_.h"
#include <assert.h>

#define LOUD 0

flint_rand_t rst;

void test1(mp_limb_t log2,const fmpz_t b,const fmpz_t z)
 {
  if( (fmpz_cmp_ui(z,0)<=0) || (log2==0) )
   return;
  int r0=fmpz_cmpabs(z,b);
  int r1=cmp_positive_log2(z,log2);
  #if LOUD
   fmpz_hex_print("2**y=",b,0);
   fmpz_hex_print("    z=",z,0);
   flint_printf("    r0/1: %d/%d\n",r0,r1);
  #endif
  if( r0>0 )
   assert( r1>0 );
  if( r0<0 )
   assert( r1<0 );
  if( r0==0 )
   assert( 0 == r1 );
 }

void
test1_ui(mp_limb_t log2,const fmpz_t b,mp_limb_t mm)
 {
  fmpz_t m; fmpz_init_set_ui(m,mm);
  test1(log2,b,m);
  fmpz_clear(m);
 }

void
test0(slong b)
 {
  fmpz_t M; fmpz_init_set_ui(M,1);
  fmpz_mul_2exp(M,M,(ulong)b);
  test1_ui(b,M,0);
  test1_ui(b,M,1);
  test1_ui(b,M,(mp_limb_t)-1);
  test1(b,M,M);                      // M
  fmpz_t n; fmpz_init_set(n,M);
  fmpz_add_ui(n,M,1);
  test1(b,M,n);                      // M+1
  fmpz_sub_ui(n,M,1);
  test1(b,M,n);                      // M-1
  slong i;
  fmpz_t m; fmpz_init(m);
  fmpz_t Mhalf; fmpz_init(Mhalf);
  fmpz_fdiv_q_2exp(Mhalf,M,1);
  for(i=100;i--;)
   {
    fmpz_randm(m,rst,n);                      // m = random(M-1)
    test1(b,M,m);
    fmpz_add(m,m,Mhalf);                     // m+M/2
    test1(b,M,m);
   }
  fmpz_mul_2exp(m,M,1);
  test1(b,M,m);                      // M<<1
  fmpz_mul_2exp(m,m,1);
  test1(b,M,m);                      // M<<2
  fmpz_mul_2exp(m,m,1);
  test1(b,M,m);                      // M<<3
  fmpz_mul(m,m,M);
  test1(b,M,m);                      // M*(M<<3)
  fmpz_clear(Mhalf);
  fmpz_clear(m);
  fmpz_clear(n);
  fmpz_clear(M);
 }

int main()
 {
  flint_randinit(rst);
  slong i;
  for(i=1;i<4*FLINT_BITS+5;i++)
   test0(i);
  flint_printf("Test passed\n");
  return 0;
 }
