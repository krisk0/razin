// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/fmpz.h>

void
fmpz_set_mpfr(fmpz_t r, const mpfr_t s, mpfr_rnd_t d)
// TODO: if t is big then there is no need to copy of t, then deallocate it
 {
  mpz_t t; mpz_init(t);
  mpfr_get_z(t, s, d);
  fmpz_set_mpz(r,t);
  mpz_clear(t);
 }
