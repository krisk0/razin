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

#include "../ulong_extras/profile_.h"

void
fmpz_mat_det_4block(fmpz_t r, const fmpz_mat_t A)
 {
  slong dim = A->r;
  if(dim < 5)
   fmpz_mat_det_cofactor(r, A);
  else
   {
    fmpz_t t; fmpz_init(t);
    // TODO: for big and fat matrice, use IML nonsingSolvLlhsMM() instead of
    //  fmpz_mat_solve_dixon()
    MARK_TIME(t0);
    fmpz_mat_det_divisor(t, A);
    DUMP_TIME("fmpz_mat_det_divisor()",t0);
    if(fmpz_is_zero(t))
     fmpz_set_ui(r,0);
    else
     {
      MARK_TIME(t1);
      fmpz_mat_det_modular_given_divisor_4block(r, A, t);
      DUMP_TIME("fmpz_mat_det_modular_given_divisor_4block()",t1);
     }
    fmpz_clear(t);
   }
 }
