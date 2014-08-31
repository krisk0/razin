// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

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
  //assert(r < p_deg_t); this assert fails
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
