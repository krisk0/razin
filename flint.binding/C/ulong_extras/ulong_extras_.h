// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef ULONG_EXTRAS__H
#define ULONG_EXTRAS__H

#include <flint/nmod_vec.h>
#include "longlong_.h"

#if !defined(SPEEDUP_NMOD_RED3)
 #define SPEEDUP_NMOD_RED3 1
#endif

typedef struct
 {
   mp_limb_t p;
   ulong k;
   mp_limb_t p_deg_k;
 } p_k_pk_t;

mp_limb_t n_pow_speedup(mp_limb_t n, ulong exp);
mp_limb_t n_mulmod_preinv_4arg(mp_limb_t a, mp_limb_t b, mp_limb_t n,
  mp_limb_t ninv);
mp_limb_t inv_mod_pk(mp_limb_t a,mp_limb_t p,ulong t,
  mp_limb_t p_deg_t,mp_limb_t p_deg_t_norm,mp_limb_t p_deg_t_inv);
mp_limb_t inv_mod_pk_3arg(mp_limb_t a,const p_k_pk_t pp,const nmod_t nn);
void invert_permutation(mp_limb_t* r,mp_limb_t* s,slong n);
mp_limb_t count_primes_in_range(mp_limb_t lo,mp_limb_t up);
void loud_primes_count(mp_limb_t lo);

static __inline__
void max_degree(p_k_pk_t* s)
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

static __inline__
void init__p_k_pk__and__nmod(p_k_pk_t* s,nmod_t* m)
// Based on s->p, initialize all other fields of s and m, 
 {
  mp_limb_t t;
  max_degree(s);
  // looks like count_leading_zeros() makes into same code as 
  //  count_leading_zeros_opt()
  count_leading_zeros_opt( t, s->p_deg_k );
  m->n = s->p_deg_k << t;
  invert_limb(m->ninv, m->n);
  #if SPEEDUP_NMOD_RED3
   t = - m->n;
   m->norm = n_mulmod_preinv_4arg( t,t, m->n,m->ninv );
  #endif
 }

static __inline__
void init_nmod_from_pp(nmod_t* m,const p_k_pk_t* s)
// Based on s, initialize m
 {
  mp_limb_t t;
  count_leading_zeros_opt( t, s->p_deg_k );
  m->n = s->p_deg_k << t;
  invert_limb(m->ninv, m->n);
  #if SPEEDUP_NMOD_RED3
   t = - m->n;
   m->norm = n_mulmod_preinv_4arg( t,t, m->n,m->ninv );
  #endif
 }

void mod_flint_style(nmod_t* m,const p_k_pk_t* s)
// change m from pk-style to regular FLINT style, to match modulo s->p
 {
  if(s->k > 1)
   {
    // no correct fields at all
    nmod_init(m,s->p);
    return;
   }
  if( s->p >= (UWORD(1)<<(FLINT_BITS-1)) )
   {
    // only m->norm is incorrect
    m->norm=0;
    return;
   }
  nmod_init(m, s->p);
 }

typedef struct
 {
  slong index;
  slong allocated_size;
  slong last_output_mod_30;
  mp_limb_t* numbers;
 }
n_primes_rev_struct;
typedef n_primes_rev_struct n_primes_rev_t[1];

int 
is_degree_of_2(mp_limb_t x)
 {
  return !( x & (x-1) );
 }

static __inline__ mp_limb_t
n_gcd_ui_positive(mp_limb_t* u,mp_limb_t* v,mp_limb_t x,mp_limb_t y,mp_limb_t n)
//return g=gcd(x,y) and set numbers u,v such that u*x+v*y=g modulo n
 {
  mp_limb_t g;
  if(x<y)
   {
    // TODO: n_xgcd is in C and GMP mpn_gcd_1 is in ASM.
    // Should I use a GMP subroutine instead of n_xgcd?
    g=n_xgcd(v,u,y,x);
    //assert( g == (*v)*y - (*u)*x );
    *u = n-( (*u) % n );
    *v %= n;
    return g;
   }
  g=n_xgcd(u,v,x,y);
  //assert( g == (*u)*x - (*v)*y );
  *v = n-( (*v) % n );
  *u %= n;
  return g;
 }

static __inline__ mp_limb_t
n_gcd_ui_positive_4arg(mp_limb_t* u,mp_limb_t* v,mp_limb_t x,mp_limb_t y)
//return g=gcd(x,y) and set numbers u,v such that u*x+v*y=g modulo T
 {
  mp_limb_t g;
  if(x<y)
   {
    // TODO: n_xgcd is in C and GMP mpn_gcd_1 is in ASM.
    // Should I use a GMP subroutine instead of n_xgcd?
    g=n_xgcd(v,u,y,x);
    //assert( g == (*v)*y - (*u)*x );
    *u = -*u;
    return g;
   }
  g=n_xgcd(u,v,x,y);
  //assert( g == (*u)*x - (*v)*y );
  *v = -*v;
  return g;
 }

#endif
