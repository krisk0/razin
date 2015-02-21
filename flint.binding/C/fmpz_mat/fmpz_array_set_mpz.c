// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

void fmpz_array_set_mpz(fmpz* tgt,mpz_t* sou,slong siz)
 {
  slong i;
  for(i=0;i<siz;i++)
   fmpz_set_mpz(tgt+i,sou[i]);
 }
