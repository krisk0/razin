// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>

typedef struct
 {
   mp_limb_t p;
   ulong k;
   mp_limb_t p_deg_k;
 } p_k_pk_t;

static __inline__
void max_degree_slow(p_k_pk_t* s)
// Initialize k and p_deg_k fields of argument, selecting max possible k
 {
  mp_limb_t p=s->p;
  // TODO: get rid of division
  mp_limb_t b=UWORD_MAX / p;
  s->k=1;
  s->p_deg_k=p;
  while( s->p_deg_k <= b )
   {
    s->p_deg_k *= p;
    ++s->k;
   }
 }

static __inline__
void max_degree_fast(p_k_pk_t* s)
// Initialize k and p_deg_k fields of argument, selecting max possible k
 {
  mp_limb_t p=s->p;
  mp_limb_t H,L,p_deg_k=p,k=1;
  while(1)
   {
    umul_ppmm( H,L, p_deg_k,p );
    if(H)
     {
      s->k=k;
      s->p_deg_k=p_deg_k;
      return;
     }
    ++k;
    p_deg_k=L;
   }
 }

void test_p(mp_limb_t p)
 {
  p_k_pk_t S0,S1;
  S0.p=p | 1;
  S1.p=p | 1;
  max_degree_slow(&S0);
  max_degree_fast(&S1);
  assert( S0.p_deg_k == S1.p_deg_k );
  assert( S0.k == S1.k );
 }

int main()
 {
  test_p(3); test_p(5);
  mp_limb_t one=1,i;
  for(i=2;i<63;i++)
   test_p(one << i);
  flint_printf("Test passed\n");
 }
