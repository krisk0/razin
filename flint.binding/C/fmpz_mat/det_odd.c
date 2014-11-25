// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>
#include <time.h>
#include <flint/flint.h>
#include "../ulong_extras/profile_.h"
#include "../fmpz_mat/fmpz_mat_.h"

#include "det_divisor_odd.c"
mp_limb_t fmpz_mat_det_divisor_odd(mpz_t r, mp_limb_t* det_mod_T, fmpz_mat_t);

void fmpz_mat_det_odd(fmpz_t r,const fmpz_mat_t a)
 {
  #if 0
   flint_printf("fmpz_mat_det_odd() A:\n");
   fmpz_mat_print_pretty(a);
   abort();
  #endif
  mp_limb_t det_mod_T;
  mpz_t m; mpz_init_set_ui(m, 1);
  MARK_TIME(t0);
  mp_limb_t hb=fmpz_mat_det_divisor_odd(m, &det_mod_T, a);
  DUMP_TIME("_det_divisor_odd()",t0);
  if( hb )
   {
    MARK_TIME(t1);
    fmpz_mat_det_modular_given_divisor_4arg(m, hb, det_mod_T, a);
    DUMP_TIME("det_modular_given_divisor_4arg()",t1);
   }
  fmpz_set_mpz(r,m);
  mpz_clear(m);
 }

#undef NDEBUG
