// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>
#include <flint/flint.h>
#include "C/ulong_extras/longlong_.h"
#include "C/ulong_extras/ulong_extras_.h"
#include "C/ulong_extras/mulmod_preinv_4arg.c"
#include "C/nmod_mat/nmod_mat_.h"

flint_rand_t st;
p_k_pk_t P;
nmod_t M;
mp_limb_t Inp[8];

mp_limb_t R()
 {
  return n_randlimb(st) % M.n;
 }

#define DOT_slow(t, s, mod)                           \
  t=n_mulmod_preinv_4arg( s[0],s[1], mod.n,mod.ninv ); \
  t=n_addmod(                                           \
    t,                                                   \
    n_mulmod_preinv_4arg( s[2],s[3], mod.n,mod.ninv ),   \
    mod.n);                                              \
  t=n_addmod(                                            \
    t,                                                   \
    n_mulmod_preinv_4arg( s[4],s[5], mod.n,mod.ninv ),   \
    mod.n);                                              \
  t=n_addmod(                                            \
    t,                                                   \
    n_mulmod_preinv_4arg( s[6],s[7], mod.n,mod.ninv ),  \
    mod.n);

#define DOT_fast(r, s, mod)                           \
  {                                            \
   VECTOR_DOT_HEAD(s[0],s[1]);                 \
   VECTOR_DOT_BODY(s[2],s[3]);                  \
   VECTOR_DOT_BODY(s[4],s[5]);                   \
   VECTOR_DOT_BODY(s[6],s[7]);                    \
   VECTOR_DOT_TAIL(r,mod.n,mod.ninv,mod.norm);    \
  }

void DOT_test()
 {
  mp_limb_t r0,r1;
  DOT_slow(r0, Inp, M);
  DOT_fast(r1, Inp, M);
  if(r0 != r1)
   {
    flint_printf("n=%lX p=%lX k=%lX\n",M.n,P.p,P.k);
    flint_printf("Inp=%lX %lX %lX %lX %lX %lX %lX %lX\n",Inp[0],Inp[1],
     Inp[2],Inp[3],Inp[4],Inp[5],Inp[6],Inp[7]);
    flint_printf("good=%wu baad=%wu\n",r0,r1);
    assert(0);
   }
 }

void test_p()
 {
  slong i,j;
  init__p_k_pk__and__nmod( &P, &M );
  #if 0
   A8C6C171C479E812*D3794110098872AF 99274D68DB68BF16*D4F5F1EA6E50110A 
   C720260C6D97F198*755B8DE06C7C50C0 DA1C83D6A9476672*B8BBF67A6D33643E 
    = 8DD926D4A87B1004
  #endif
  Inp[0] = 0xA8C6C171C479E812 % M.n;
  Inp[1] = 0xD3794110098872AF % M.n;
  Inp[2] = 0x99274D68DB68BF16 % M.n;
  Inp[3] = 0xD4F5F1EA6E50110A % M.n;
  Inp[4] = 0xC720260C6D97F198 % M.n;
  Inp[5] = 0x755B8DE06C7C50C0 % M.n;
  Inp[6] = 0xDA1C83D6A9476672 % M.n;
  Inp[7] = 0xB8BBF67A6D33643E % M.n;
  DOT_test();
  for(i=100;i--;)
   {
    for(j=8;j--;)
     Inp[j]=R();
    DOT_test();
   }
 }

int main()
 {
  flint_randinit(st);
  P.p=3;

  while(1)
   {  
    test_p();
    P.p=n_nextprime(P.p,0);
    if(P.p>500)
     break;
   }
  
  flint_printf("Test passed\n");
 }
