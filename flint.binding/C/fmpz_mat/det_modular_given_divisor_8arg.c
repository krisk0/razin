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
#undef NDEBUG
#include <assert.h>

static __inline__ void
decrease_bound_fmpz(mpfr_t b,mpfr_prec_t pr,mpz_t d)
 {
  mpfr_t dF; mpfr_init2(dF, mpz_sizeinbase(d,2));
  mpfr_t log2_d; mpfr_init2(log2_d,pr);

  mpfr_set_z(dF, d, MPFR_RNDZ);
  mpfr_log2(log2_d, dF, MPFR_RNDZ);
  mpfr_sub(b, b, log2_d, MPFR_RNDU);

  mpfr_clear(log2_d);
  mpfr_clear(dF);
 }

static __inline__ void
decrease_bound_ui(mpfr_t b,mpfr_prec_t pr,mp_limb_t d)
 {
  mpfr_t dF; mpfr_init2(dF, FLINT_BITS);
  mpfr_t log2_d; mpfr_init2(log2_d,pr);

  mpfr_set_uj(dF, d, MPFR_RNDZ);
  mpfr_log2(log2_d, dF, MPFR_RNDZ);
  mpfr_sub(b, b, log2_d, MPFR_RNDU);

  mpfr_clear(log2_d);
  mpfr_clear(dF);
 }

static __inline__ int
comp_bound_ui(mpfr_t b,mp_limb_t q)
// 1 iff log2(q) rounded down <= b
 {
  mpfr_t log2; mpfr_init(log2);
  mpfr_t qF; mpfr_init2(qF, FLINT_BITS);
  mpfr_set_uj(qF, q, MPFR_RNDZ);
  mpfr_log2(log2, qF, MPFR_RNDZ);
  int r=(mpfr_cmp(log2,b)<0);
  mpfr_clear(qF);
  mpfr_clear(log2);
  return r;
 }

static __inline__ mp_limb_t
invert_det_divisor_modulo_pk(mpz_t dd,p_k_pk_t const* pp,nmod_t const* mod)
// take dd modulo p_deg_k, then invert it
 {
  mp_limb_t r,m=pp->p_deg_k
  r=mpz_fdiv_ui(divisor, m);
  if(pp->k==1)
   return n_invmod(r,m);
  else
   return inv_mod_pk_3arg(r,pp[0],mod[0]);
 }

static __inline__ void
mpz_fmpz_mul_2arg(mpz_t z,fmpz_t x)
 {
  register slong xx=(slong)(*x);
  if(!COEFF_IS_MPZ(xx))
   {                              // x is small
    if(xx != WORD(1))
     mpz_mul_ui(z,z,(mp_limb_t)xx);
   }
  else
   mpz_mul(z,z,COEFF_TO_PTR(xx)); // x is big
 }

void 
fmpz_mat_det_modular_given_divisor_8arg(mpz_t det, nmod_mat_t Amod,
  mpfr_t hadamard_log2, mpfr_prec_t pr, p_k_pk_t* pp, n_primes_rev_t iT,
  mp_limb_t xmod, const fmpz_mat_t A)
/*
act like fmpz_mat_det_modular_given_divisor_4block(), but
 * decrease primes using iT, rather than increase them
 * don't count H.B., use hadamard_log2 which is 1+log2(H.B.)
 * sum logarithms instead of multiplying together found primes
 * re-use found prime pp->p and xmod which is determinant of A modulo pp->p

hadamard_log2 on entry is upper bound on log2(2*H.B)
              decreased to a non-positive value on exit

iT on entry just found prime pp->p
   possibly gets shifted on exit
*/
 {
  // loop bound = 2*H.B / known det divisor
  decrease_bound_fmpz(hadamard_log2,pr,det);

  // re-use known det A modulo pp->p_deg_k
  mp_limb_t divisor_inv=invert_det_divisor_modulo_pk(det,pp,&Amod->mod);
  xmod=n_mulmod2_preinv(xmod,divisor_inv, Amod->mod.n,Amod->mod.ninv);
  fmpz_t xnew,x;
  fmpz_init(xnew); 
  fmpz_init(x); 
  fmpz_t prod; mpz_init_set_ui( prod, UWORD(1) );

  fmpz_CRT_ui(xnew, x, prod, xmod, pp->p_deg_k, 1);
  fmpz_set_ui(prod, pp->p_deg_k);
  fmpz_set(x, xnew);

  // for orthogonal matrice the bound might be reached at this point.
  //  Attempt to skip main loop
  if(comp_bound_ui(hadamard_log2,pp->p_deg_k))
   {
    mp_limb_t* scratch=flint_malloc( 4*(A->r-4)*sizeof(mp_limb_t) );
    mp_limb_t bound=mpfr_get_uj(hadamard_log2,MPFR_RNDU);

    for(;;)
     {
      divisor_inv=choose_prime_and_degree( &pp, &Amod->mod, it, det );
      mpz_mat_get_nmod_mat(Amod, A);
      xmod=nmod_mat_det_mod_pk_4block(Amod,pp,scratch);
      xmod=n_mulmod2_preinv(xmod,divisor_inv, Amod->mod.n,Amod->mod.ninv);
      fmpz_CRT_ui(xnew, x, prod, xmod, pp->p_deg_k, 1);
      fmpz_mul_ui(prod, prod, pp->p_deg_k);
      if(cmp_positive_log2(prod,bound) >= 0)
       break;
      fmpz_set(x, xnew);
     }
  
    flint_free(scratch);
   }

  fmpz_clear(prod);
  mpz_fmpz_mul_2arg(det,xnew);
  fmpz_clear(prod);
  fmpz_clear(x);
  fmpz_clear(xnew);
 }

#undef NDEBUG
