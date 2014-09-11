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

#endif
