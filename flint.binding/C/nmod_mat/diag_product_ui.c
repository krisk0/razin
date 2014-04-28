// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// A FLINT developer called trace smth different than my teacher of algebra.
//  So call it diag_product

#include <flint/nmod_mat.h>

mp_limb_t nmod_mat_diag_product_ZZ_ui(const nmod_mat_t m)
// counts product of diagonal entries (treated as integers) modulo 2**64
// m: square, dimension>0
 {
  slong cc=m->c-1;
  mp_limb_t** rows=m->rows;
  mp_limb_t r=rows[cc][cc];
  for( ; cc--; )
   r *= rows[cc][cc];
  return r;
 }
