// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef FMPZ_MAT__H
#define FMPZ_MAT__H

#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include "../fmpz/fmpz_.h"
#include "../nmod_mat/nmod_mat_.h"
#include "../mpz_square_mat/mpz_square_mat_.h"

void fmpz_mat_det_modular_given_divisor_8arg(mpz_t det, nmod_mat_t Amod,
  mpfr_t hadamard_log2, mpfr_prec_t pr, p_k_pk_t* pp, n_primes_rev_t it,
  mp_limb_t xmod, const fmpz_mat_t A);
void fmpz_mat_det_suspected_zero(mpz_t r,const fmpz_mat_t A,const mpz_t W);
void fmpz_triU_inverse_smallDet(fmpz_mat_t r, fmpz_t d, const fmpz_mat_t s);
void fmpz_mat_scalar_divexact_ui_2arg(fmpz_mat_t R,mp_limb_t d);
void fmpz_mat_det_odd(fmpz_t r,const fmpz_mat_t a);
slong hadamard_2arg(mpfr_t b,const fmpz_mat_t m);
mp_limb_t cramer_rule(const mpfr_t den_bound, 
 mpz_square_mat_t A, mpfr_prec_t pr, slong k);
int fmpz_mat_hermitian_decomposition_2(fmpz_mat_t b,fmpz_t r,
 const fmpz_mat_t m);

#endif
