// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>
#include <flint/flint.h>
#include "../fmpz_mat/fmpz_mat_.h"

void fmpz_mat_det_odd(fmpz_t r,const fmpz_mat_t A)
 {
  flint_printf("fmpz_mat_det_odd() A:\n");
  fmpz_mat_print_pretty(A);
  abort();
 }

#undef NDEBUG
