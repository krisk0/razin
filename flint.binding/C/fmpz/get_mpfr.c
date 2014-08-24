// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/fmpz.h>
#include "fmpz_.h"

void
fmpz_get_mpfr_3arg(mpfr_t r, fmpz s, mpfr_rnd_t d)
 {
  fmpz_get_mpfr_slave(r,s,d);
 }

void
fmpz_get_mpfr(mpfr_t r, const fmpz_t s, mpfr_rnd_t d)
 {
  fmpz_get_mpfr_macro(r,s[0],d);
 }
