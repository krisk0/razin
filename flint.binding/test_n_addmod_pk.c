// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>
#include <flint/flint.h>
#include "flint/ulong_extras.h"
#include "C/ulong_extras/longlong_.h"

#if 0
This program tests n_addmod_pk() macro against FLINT n_addmod()
#endif

flint_rand_t st;

void 
test_triple(mp_limb_t a,mp_limb_t b,mp_limb_t n)
 {
  assert(a<n);
  assert(b<n);
  mp_limb_t r=a;
  n_addmod_pk(r,b,n);
  mp_limb_t x=n_addmod(a,b,n);
  assert(x<n);
  assert(r<n);
  if( x == r )
   return;
  flint_printf("\n(%wu + %wu) mod %wu = %wu != %wu\n",a,b,n,x,r);
  assert(0);
 }

void 
test_prime(mp_limb_t p)
 {
  p |= UWORD( 0x8000000000000000 );
  test_triple( n_randbits(st,64) % p, n_randbits(st,64) % p, p );
 }
 
int main()
 {
  flint_randinit(st);
  test_prime(3);
  test_prime(0xFFFFFFFFFFFFFFC5);
  slong i;
  for(i=100;i--;)
   test_prime( n_randbits(st,64) );
  flint_printf("Test passed\n");
 }
