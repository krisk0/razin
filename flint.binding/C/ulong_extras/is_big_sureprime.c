// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

/*
FLINT n_is_probabprime() is not very good for numbers >= 10**16

The subroutine below should be better. It should test primality without errors,
 if there is no mistake on http://miller-rabin.appspot.com concerning base 
2 325 9375 28178 450775 9780504 1795265022
*/

#include "is_strong_probabprime2_preinv.c"

static __inline__ int 
n_is_big_sureprime(mp_limb_t n)
/*
return 1 iff MR test passes for a=2 325 9375 28178 450775 9780504 1795265022
 
if n<10**16, then n_is_probabprime() should be faster
*/
 {
  mp_limb_t d = n-1;
  unsigned int norm;
  mp_limb_t ninv;
  count_trailing_zeros(norm, d);
  d >>= norm;
  ninv = n_preinvert_limb(n);
  #define MR n_is_strong_probabprime2_preinv_speedup
  return MR(n, ninv, UWORD(2), d) && 
         MR(n, ninv, UWORD(325), d) &&
         MR(n, ninv, UWORD(9375), d) &&
         MR(n, ninv, UWORD(28178), d) &&
         MR(n, ninv, UWORD(450775), d) &&
         MR(n, ninv, UWORD(9780504), d) &&
         MR(n, ninv, UWORD(1795265022), d);
  #undef MR
 }
