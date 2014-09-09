// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

#define RAT_REC_TAKES_D_SERIOUSLY 1

void det_divisor_rational_reconstruction(mpz_t d,mpz_ptr x,mpz_t M,mp_limb_t p,
  slong n,mp_limb_t log2_N, mp_limb_t log2_D);
