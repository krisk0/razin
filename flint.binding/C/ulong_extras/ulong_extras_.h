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

#endif
