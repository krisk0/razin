// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz_mat.h>
#include <mpfr.h>
#include "../ulong_extras/ulong_extras_.h"
#undef NDEBUG
#include <assert.h>
void nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim);
mp_limb_t nmod_mat_det_mod_pk_4block(nmod_mat_t M,const p_k_pk_t pp,mp_limb_t* scrtch);

static __inline__ mp_limb_t
choose_prime_degree(p_k_pk_t* pp,nmod_t* mod,const fmpz_t divisor,
  n_primes_rev_t it)
 {
  assert(0);
 }

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
  assert(0);

  slong dim = A->r;
  mp_limb_t* scratch=flint_malloc( 4*(dim-4)*sizeof(mp_limb_t) );
  #if !SPEEDUP_NMOD_RED3
   Amod->mod.norm=0;
  #endif
  
  flint_free(scratch);
 }

#undef NDEBUG
