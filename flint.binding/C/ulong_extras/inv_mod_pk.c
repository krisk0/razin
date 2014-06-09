// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2009, 2010, 2012 William Hart

******************************************************************************/
// Copyright Денис Крыськов 2014

// This file contains a pieve of FLINT mulmod_preinv.c owned by W.Hart

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

static mp_limb_t
n_mulmod_preinv_4arg(mp_limb_t a, mp_limb_t b, mp_limb_t n, 
                                          mp_limb_t ninv)
//castrated FLINT subroutine n_mulmod_preinv(): norm parameter removed
{
    mp_limb_t q0, q1, r, p_hi, p_lo;

    /* multiply */
    umul_ppmm(p_hi, p_lo, a, b);
    
    /* reduce mod n */
    {
        umul_ppmm(q1, q0, ninv, p_hi);
        add_ssaaaa(q1, q0, q1, q0, p_hi, p_lo);

        r = (p_lo - (q1 + 1) * n);

        if (r >= q0)
            r += n;

        return (r < n ? r : r - n);
    }
}

mp_limb_t
inv_mod_pk(mp_limb_t a,mp_limb_t p,ulong t,
  mp_limb_t p_deg_t,mp_limb_t p_deg_t_norm,mp_limb_t p_deg_t_inv)
/*
return inverse of a modulo p**t where p is odd prime and p**t < 2**64, p does
 not divide a
*/
 {
  mp_limb_t r=n_invmod(a,p),q;
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
  return r;
 }
