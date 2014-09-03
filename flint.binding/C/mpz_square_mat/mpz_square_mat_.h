// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef MPZ_SQUARE_MAT__H
#define MPZ_SQUARE_MAT__H

#include <flint/fmpz_mat.h>

typedef struct
 {
  mpz_ptr entries;
  slong r;
  mpz_ptr* rows;
  slong* size_in_limbs;
 }
mpz_square_mat_struct;

typedef mpz_square_mat_struct mpz_square_mat_t[1];

void mpz_square_mat_init_set(mpz_square_mat_t R,const fmpz_mat_t S);
void mpz_square_mat_clear();

#endif
