// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/
// Copyright Денис Крыськов 2014

// This file contains pieces of FLINT det*.c owned by F.Johansson

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
void fmpz_mat_det_modular_given_divisor_4block(fmpz_t r,const fmpz_mat_t A,
  fmpz_t divisor);

void
fmpz_mat_det_4block(fmpz_t r, const fmpz_mat_t A)
 {
  slong dim = A->r;
  if(dim < 5)
   fmpz_mat_det_cofactor(r, A);
  else
   {
    fmpz_t t; fmpz_init(t);
    fmpz_mat_det_divisor(t, A);
    fmpz_mat_det_modular_given_divisor_4block(r, A, t);
    fmpz_clear(t);
   }
 }
