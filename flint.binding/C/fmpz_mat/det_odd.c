// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>
#include <time.h>
#include <flint/flint.h>
#include "../fmpz_mat/fmpz_mat_.h"

//#include "det_divisor_odd.c"

void fmpz_mat_det_odd(fmpz_t r,const fmpz_mat_t a)
 {
  flint_printf("fmpz_mat_det_odd() A:\n");
  fmpz_mat_print_pretty(a);
  abort();
  mp_limb_t det_mod_T;
  mpz_t m; mpz_init(m);
  #if 0
  mp_limb_t hb=fmpz_mat_det_divisor_odd(m, &det_mod_T, a);
  if( hb )
   fmpz_mat_det_modular_given_divisor_4arg(m, hb, det_mod_T, a);
  #endif
  fmpz_set_mpz(r,m);
  mpz_clear(m);
 }

#undef NDEBUG
