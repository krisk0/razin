// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

void
det_divisor_rational_reconstruction(mpz_t d,mpz_ptr x,mpz_t M,mp_limb_t p,
  slong n,mp_limb_t log2_N, mp_limb_t log2_D)
/*
M=modulo=p**k
n: length of x
log2_N: log2(upper bound on numerators)
log2_D: log2(upper bound on denominators)
*/
 {
  mpz_set_ui(d,1);
  flint_printf("det_divisor_rational_reconstruction not imple\n");
  abort();
 }
