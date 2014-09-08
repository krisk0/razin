// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef FMPZ_MAT__H
#define FMPZ_MAT__H

#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include "../nmod_mat/nmod_mat_.h"

void fmpz_mat_det_modular_given_divisor_8arg(mpz_t det, nmod_mat_t Amod,
  mpfr_t hadamard_log2, mpfr_prec_t pr, p_k_pk_t* pp, n_primes_rev_t it,
  mp_limb_t xmod, const fmpz_mat_t A);

#endif
