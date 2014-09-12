// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/fmpz_mat.h>

void
fmpz_mat_scalar_divexact_ui_2arg(fmpz_mat_t R,mp_limb_t d)
 {
  fmpz* q=R->entries;
  slong s=R->r * R->c;
  for(;s--;q++)
   fmpz_divexact_ui( q, q, d );
 }
