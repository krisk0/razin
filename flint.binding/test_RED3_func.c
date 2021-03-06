// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>
#include <string.h>                             // for memcpy
#include <flint/flint.h>
#include "C/ulong_extras/mulmod_preinv_4arg.c"
#include <flint/nmod_vec.h>
#include "C/ulong_extras/longlong_.h"
#include "C/ulong_extras/ulong_extras_.h"

#if 0
This program tests function NMOD_RED3_pk_func() and macro NMOD_RED3_pk against
 macro NMOD_RED3
#endif

flint_rand_t st;

void 
test_prime(mp_limb_t p)
 {
  p_k_pk_t P;
  P.p=p;
  flint_printf("p=%wu\n",p);
  nmod_t M0,M1;
  init__p_k_pk__and__nmod( &P, &M0 );
  memcpy( &M1, &M0, sizeof(M1) );
  assert(M0.norm);
  M0.norm=0;
  slong i,j;
  mp_limb_t S0[3],S1[3],S2[3],S3[3];
  mp_limb_t r;
  for(i=33554432;i--;)
   {
    for(j=3;j--;)
     S3[j]=S2[j]=S1[j]=S0[j]=n_randlimb(st);
    S3[2] %= M1.n;
    NMOD_RED( S0[2], S0[2], M0 );
    NMOD_RED3( S0[0], S0[2], S0[1], S0[0], M0 );
    r = NMOD_RED3_pk_func( S1[2],S1[1],S1[0], M1 );
    assert( r==S0[0] );
    NMOD_RED3_pk( S2[2],S2[1],S2[0], M1.n,M1.ninv,M1.norm );
    assert( r==S2[0] );
    NMOD_RED3_pk_5arg( S3[2],S3[1],S3[0], M1.n,M1.ninv );
    assert( r==S3[0] );
   }
 }
 
int main()
 {
  flint_randinit(st);
  test_prime( 4294967291 );
  test_prime(3);
  //test_prime(5);     5 will be tested for i=2
  test_prime(0xFFFFFFFFFFFFFFC5);
  slong i;
  mp_limb_t one=1;
  for(i=2;i<64;i++)
   {
    test_prime( n_nextprime( one << i,         0) );
    test_prime( n_nextprime( n_randbits(st,i), 0) );
   }
  test_prime( n_nextprime( n_randbits(st,64), 0) );
  flint_printf("Test passed\n");
 }
