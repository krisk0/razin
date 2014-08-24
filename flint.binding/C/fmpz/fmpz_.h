// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef FMPZ__H
#define FMPZ__H

#define PRINTF

// s should be fmpz, s evaluated multiple times
#define fmpz_get_mpfr_slave(r, s, rnd)    \
 if( COEFF_IS_MPZ(s) )                      \
   (void)mpfr_set_z(r, COEFF_TO_PTR(s), rnd); \
 else                                         \
   mpfr_set_si(r, s, rnd);                   

// _s should be fmpz, _s evaluated once
#define fmpz_get_mpfr_macro(r, _s, rnd) \
 {                                     \
  register fmpz _sR=_s;               \
  fmpz_get_mpfr_slave(r,_sR,rnd);    \
 }

#endif
