// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef ULONG_EXTRAS__H
#define ULONG_EXTRAS__H

#include <flint/nmod_vec.h>

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
void init__p_k_pk__and__nmod(p_k_pk_t* s,nmod_t* m)
// Initialize all fields of s and m, based on s->p
 {
  mp_limb_t t;
  max_degree(s);
  count_leading_zeros( t, s->p_deg_k );
  m->n = s->p_deg_k << t;
  m->ninv = n_preinvert_limb(m->n);
 }

#endif
