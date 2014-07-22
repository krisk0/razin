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
  ulong k=1;
  ulong t_div_2=t>>1;
  mp_limb_t p_deg_k=p,p_deg_k_plus_m;
  while(k <= t_div_2)
   {
    p_deg_k_plus_m=p_deg_k*p_deg_k;
    if( p_deg_k_plus_m < UWORD(1)<<(FLINT_BITS/2) )
     // can't overflow because a,r < p**(2*k) < 2**32
     q=(a%p_deg_k_plus_m)*r % p_deg_k_plus_m - 1;
    else
     {
      // TODO: is it faster to use n_mulmod2_preinv() instead?
      q=n_mulmod_preinv_4arg(a,r,p_deg_t_norm,p_deg_t_inv); 
      q = (q % p_deg_k_plus_m)-1;
     }
    q /= p_deg_k; // should divide exactly
    q = q*r % p_deg_k; // no overflow because q,r < p**k < 2**32
    r += p_deg_k * (p_deg_k-q);
    assert(r < p_deg_t);
    k <<= 1; p_deg_k = p_deg_k_plus_m;
   }
  if(k < t)
   {
    // TODO: is it faster to use n_mulmod2_preinv() instead?
    q=(n_mulmod_preinv_4arg(a,r,p_deg_t_norm,p_deg_t_inv)-1) % p_deg_t; //
    q /= p_deg_k; // should divide exactly
    p_deg_k_plus_m = p_deg_t / p_deg_k;
    // m < t/2 => p**m < 2**32, no overflow in next line
    q = q * (r % p_deg_k_plus_m) % p_deg_k_plus_m;
    r += p_deg_k*(p_deg_k_plus_m-q);
   }
  return r % p_deg_t;

  #undef p
  #undef t 
  #undef p_deg_t
  #undef p_deg_t_norm 
  #undef p_deg_t_inv 
 }
