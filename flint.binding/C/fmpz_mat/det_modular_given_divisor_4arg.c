// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz_mat.h>
#include <mpfr.h>
#include "../ulong_extras/ulong_extras_.h"
#include "../fmpz/fmpz_.h"
#include "../nmod_mat/nmod_mat_.h"
#include "../tmod_mat/tmod_mat_.h"
#undef NDEBUG
#include <assert.h>

#define LOUD_DET_BOUND 1

//mpz_get_ui (const mpz t op)

static __inline__ mp_limb_t
div_modulo_T(mp_limb_t y,mp_limb_t a)
// y / a modulo 2**FLINT_BITS
 {
  if( y == a )        // this should happen often
   return UWORD(1);
  a=t_invmod(a);
  return y*a;  // hope compiler multiplies unsigned modulo 2**64 correctly
 }

static __inline__ void
_CRT_data_init_20140912(fmpz_t x0,fmpz_t x1,fmpz_t prod,mp_limb_t det_mod_T,
  const mpz_t dd)
 {
  // prod = T
  fmpz_init_set_ui(prod,UWORD(1)); fmpz_mul_2exp(prod,prod,FLINT_BITS);
  // x0 = (det_mod_T / dd) symmmod T
  mp_limb_t t=div_modulo_T(det_mod_T,flint_mpz_get_ui(dd)); // dd is positive
  fmpz_init_set_ui(x0,t);
  if( t & (UWORD(1)<<(FLINT_BITS-1)) )
   fmpz_sub(x0,x0,prod); // symmetric range CRT, must not exceed M/2
  fmpz_init_set(x1,x0);
 }

static __inline__ mp_limb_t
_choose_prime_20140912(p_k_pk_t* pp, nmod_t* mod, n_primes_rev_t it,
  const mpz_t divisor)
 {
  mp_limb_t r,t=1,r_mod_p;
  for(;;)
   {
    if( 0==pp->p )
     pp->p=n_primes_rev_init(it, 0);
    else
     {
      if( (pp->p=n_primes_rev_next(it)) == 1 )
       {
        flint_printf("Exception (fmpz_mat_det_modular_given_divisor_4arg): "
                       "Prime set exhausted\n");
        abort();
       }
     }
    if( pp->p >= (UWORD(1)<<(FLINT_BITS/2)) )
     {
      pp->k=1;
      r=r_mod_p=mpz_fdiv_ui( divisor, pp->p_deg_k=pp->p );
      t=0;
     }
    else
     {
      max_degree( pp );
      r=mpz_fdiv_ui( divisor, pp->p_deg_k );
      r_mod_p=r % pp->p;
     }
    if(r_mod_p)
     {
      if(t)
       count_leading_zeros( t, pp->p_deg_k );
      mod->n = pp->p_deg_k << t;
      invert_limb(mod->ninv, mod->n);
      #if SPEEDUP_NMOD_RED3
       t = - mod->n;
       mod->norm = n_mulmod_preinv_4arg( t,t, mod->n,mod->ninv );
      #else
       mod->norm = 0;
      #endif
      return inv_mod_pk_4arg(r,r_mod_p,pp[0],mod[0]);
     }
   }
 }

__inline__ static void
_20140912_no_CRT(mpz_t r,mp_limb_t det_mod_t)
 {
  mp_limb_t high_bit=(UWORD(1)<<(FLINT_BITS-1));
  if(mpz_cmp_ui(r,1))
   {
    // r *= symm_abs(det_mod_t / r mod t, t)
    det_mod_t=div_modulo_T(det_mod_t, flint_mpz_get_ui(r)); // r is positive
    if( det_mod_t & high_bit )
     {
      flint_mpz_mul_ui(r,r,-det_mod_t);
      mpz_neg(r,r);
     }
    else
     flint_mpz_mul_ui(r,r,det_mod_t);
   }
  else
   {
    // r = symm_abs(det_mod_T,t)
    if( det_mod_t & high_bit )
     {
      flint_mpz_set_ui(r, -det_mod_t);
      mpz_neg(r,r);
     }
    else
     flint_mpz_set_ui(r, det_mod_t);
   }
 }

void
fmpz_mat_det_modular_given_divisor_4arg(mpz_t r, mp_limb_t hb,
  mp_limb_t det_mod_T, fmpz_mat_t A)
 {
  if(hb <= FLINT_BITS)
   {// don't need Chinese Remaider or det modulo prime
    _20140912_no_CRT(r,det_mod_T);
    return;
   }
  {
   // looks like small primes are better than large primes
   // TODO: use small primes produces with n_nextprime() instead of big primes
   #if LOUD_DET_BOUND
    gmp_printf("fmpz_mat_det_modular_given_divisor_4arg(): bound=%Md\n",hb);
    slong primes_used=0;
   #endif
   p_k_pk_t pp; pp.p=0;
   n_primes_rev_t iT;
 
   fmpz_t x,xnew,prod;
   _CRT_data_init_20140912(x,xnew,prod,det_mod_T,r);
   
   mp_limb_t n=A->r;
   mp_limb_t* scratch=flint_malloc( 4*(n-4)*sizeof(mp_limb_t) );
   nmod_mat_t Amod;
   nmod_mat_init_square_2arg(Amod,n);
   
   for(;;)
    {
     mp_limb_t divisor_inv=_choose_prime_20140912(&pp, &Amod->mod, iT, r);
     fmpz_mat_get_nmod_mat(Amod, A);
     mp_limb_t xmod=nmod_mat_det_mod_pk_4block(Amod,pp,scratch);
     xmod=n_mulmod2_preinv(xmod,divisor_inv, Amod->mod.n,Amod->mod.ninv);
     fmpz_CRT_ui(xnew, x,prod, xmod,pp.p_deg_k, 1);
     fmpz_mul_ui(prod, prod, pp.p_deg_k);
     #if LOUD_DET_BOUND
      primes_used++;
     #endif
     if(cmp_positive_log2(prod,hb) >= 0)
      break;
     fmpz_set(x, xnew);
    }
 
   n_primes_rev_clear(iT);
   nmod_mat_clear(Amod);
   free(scratch);
   mpz_fmpz_mul_det_2arg(r,xnew);
   fmpz_clear(xnew);
   fmpz_clear(x);
   fmpz_clear(prod);
   #if LOUD_DET_BOUND
    flint_printf("fmpz_mat_det_modular_given_divisor_4arg() primes used: %d\n",
     primes_used);
   #endif
  }
 }

#undef NDEBUG
