// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include "ulong_extras_.h"

mp_limb_t
inv_mod_pk(mp_limb_t a,mp_limb_t p,ulong t,
  mp_limb_t p_deg_t,mp_limb_t p_deg_t_norm,mp_limb_t p_deg_t_inv)
/*
return inverse of a modulo p**t where p is odd prime and p**t < 2**64, p does
 not divide a
*/
 {
  mp_limb_t r=n_invmod(a,p),q;
  #include "det_mod_pk.inc"
 }

mp_limb_t
inv_mod_pk_3arg(mp_limb_t a,const p_k_pk_t pp,const nmod_t nn)
 {
  #define p pp.p
  #define t pp.k
  #define p_deg_t pp.p_deg_k
  #define p_deg_t_norm nn.n
  #define p_deg_t_inv nn.ninv

  mp_limb_t r=n_invmod(a,p),q;
  #include "det_mod_pk.inc"

  #undef p
  #undef t 
  #undef p_deg_t
  #undef p_deg_t_norm 
  #undef p_deg_t_inv 
 }

mp_limb_t
inv_mod_pk_4arg(mp_limb_t a,mp_limb_t a_mod_p, 
  const p_k_pk_t pp,const nmod_t nn)
 {
  #define p pp.p
  #define t pp.k
  #define p_deg_t pp.p_deg_k
  #define p_deg_t_norm nn.n
  #define p_deg_t_inv nn.ninv

  mp_limb_t r=n_invmod(a_mod_p,p),q;
  #include "det_mod_pk.inc"

  #undef p
  #undef t 
  #undef p_deg_t
  #undef p_deg_t_norm 
  #undef p_deg_t_inv 
 }
