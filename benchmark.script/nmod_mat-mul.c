// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

/*
This program tests/benchmarks RAZIN subroutine nmod_mat_mul_pk_classical() 
 against FLINT nmod_mat_mul()
*/

#include <time.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#define SPEEDUP_NMOD_RED3 1
#include <flint/flint.h>
#include "../flint.binding/C/ulong_extras/longlong_.h"
#include "../flint.binding/C/ulong_extras/ulong_extras_.h"
#include "../flint.binding/C/nmod_mat/init_3arg.c"
#include "../flint.binding/C/nmod_mat/randfill.c"
#include "../flint.binding/C/ulong_extras/mulmod_preinv_4arg.c"
#include "../flint.binding/C/nmod_mat/mul_pk_classical.c"

#define DIM_EQUAL(alpha,betta) ( (alpha->r==betta->r) && (alpha->c==betta->c) )

// multiply matrix 2*(v+1) times, put result into R1
#define CHAIN_MULT(v, f, A, B, R0, R1)    \
  update_me=R1->rows[0];                    \
  if( DIM_EQUAL(A,R1) )                      \
   {                                          \
    f( R0,  A, B );                             \
    f( R1, R0, B );                              \
    update_me[0]=1;                               \
    for(i=v;i--;)                                  \
     {                                              \
      f( R0, R1, B );                                \
      f( R1, R0, B );                                 \
      update_me[1]=1;                                  \
     }                                                  \
   }                                                     \
  else                                                    \
   {                                                       \
    if( DIM_EQUAL(B,R1) )                                   \
     {                                                       \
      f( R0, A,  B );                                         \
      f( R1, A, R0 );                                          \
      update_me[0]=1;                                           \
      for(i=v;i--;)                                              \
       {                                                          \
        f( R0, A, R1 );                                            \
        f( R1, A, R0 );                                             \
        update_me[1]=1;                                              \
       }                                                              \
     }                                                                 \
    else                                                                \
     {                                                                   \
      flint_printf("Dimension not suitable for chain multiplication\n"); \
      exit(1);                                                          \
     }                                                               \
   }                                                               

flint_rand_t st;

void
init_matrix( 
  nmod_mat_t a, nmod_mat_t b, 
  nmod_mat_t c0, nmod_mat_t c1, nmod_mat_t c2,
  slong dim0, slong dim1, slong dim2, const nmod_t M)
 {
  nmod_mat_init_3arg( a, dim0, dim1 ); memcpy( &a->mod, &M, sizeof(nmod_t) );
  nmod_mat_init_3arg( b, dim1, dim2 ); memcpy( &b->mod, &M, sizeof(nmod_t) );
  nmod_mat_randfill( a, st );
  nmod_mat_randfill( b, st );
  nmod_mat_init_3arg( c0, dim0, dim2 );  
  nmod_mat_init_3arg( c1, dim0, dim2 );  
  nmod_mat_init_3arg( c2, dim0, dim2 );  
 }

void 
reduce_modulo(nmod_mat_t t, nmod_mat_t s, const nmod_t M)
 {
  nmod_mat_init_3arg(t, s->r, s->c ); memcpy( &t->mod, &M, sizeof(nmod_t) );
  slong siz=s->r*s->c,i;
  mp_limb_t* tE=t->entries;
  mp_limb_t* sE=s->entries;
  mp_limb_t n=M.n;
  for(i=siz;i--;)
   tE[i] = sE[i] % n;
 }

void 
equality_check( nmod_mat_t a, nmod_mat_t b )
 {
  mp_limb_t n=b->mod.n;
  assert( n <= a->mod.n );
  slong siz=b->r*b->c,i;
  mp_limb_t* aE=a->entries;
  mp_limb_t* bE=b->entries;
  for(i=siz;i--;)
   assert( aE[i] % n == bE[i] );
 }

void
nmod_mat_multiplication_test_benchmark(
  slong dim0,slong dim1,slong dim2,
  mp_limb_t p,
  slong experiments)
 {
  assert(experiments >= 1);
  assert(dim2 > 1);
  clock_t t0,t1,t2,t3;
  slong i;
  mp_limb_t* update_me;
  nmod_mat_t A0,B0,C0,A1,B1,C1,C2;
  p_k_pk_t P; P.p=p;
  nmod_t M0,M1;
  flint_randinit(st);
  init__p_k_pk__and__nmod( &P, &M0 );
  nmod_init( &M1, P.p_deg_k );
  init_matrix(A0,B0, C0,C1,C2, dim0,dim1,dim2, M0);
  reduce_modulo( A1, A0, M1);
  reduce_modulo( B1, B0, M1);
  --experiments;
  memcpy( &C2->mod, &M1, sizeof(nmod_t) );
  memcpy( &C1->mod, &M1, sizeof(nmod_t) );
  t0=clock();
  CHAIN_MULT(experiments, nmod_mat_mul,    A1, B1, C2, C1);
  t1=clock();
  memcpy( &C2->mod, &M0, sizeof(nmod_t) );
  memcpy( &C0->mod, &M0, sizeof(nmod_t) );
  t2=clock();
  CHAIN_MULT(experiments, nmod_mat_mul_pk_classical, A0, B0, C2, C0);
  t3=clock();
  equality_check( C0, C1 );
  flint_printf("2**%w * %wu**%wu | %w %w %w: ",M1.norm,p,P.k,dim0,dim1,dim2);
  // ************************ Warning ************************
                                          // line below works on amd64/Linux,
  flint_printf("%w %w\n", t1-t0, t3-t2); //   because clock_t is slong
  // If clock_t is not integer, replace format string above with smth else
 }

#undef DIM_EQUAL

int main()
 {
  #define VOL 1<<18
  #define BIG_P 0xFFFFFFFFFFFFFFc5
  #define BENCHMARK nmod_mat_multiplication_test_benchmark

  flint_randinit(st);

  BENCHMARK( 4,4,100, 3, VOL );
  BENCHMARK( 100,4,4, 3, VOL );

  BENCHMARK( 4,4,100, 5, VOL );
  BENCHMARK( 100,4,4, 5, VOL );
  
  BENCHMARK( 4,4,100, BIG_P, VOL );
  BENCHMARK( 100,4,4, BIG_P, VOL );
 }

#if 0
result for my core-i5:
2**0 * 3**40 | 4 4 100: 4810000 1910000
2**0 * 3**40 | 100 4 4: 4870000 1990000
2**1 * 5**27 | 4 4 100: 3790000 2330000
2**1 * 5**27 | 100 4 4: 3780000 2200000
2**0 * 18446744073709551557**1 | 4 4 100: 4870000 2620000
2**0 * 18446744073709551557**1 | 100 4 4: 4920000 2580000
#endif
