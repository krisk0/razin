// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <gmp.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include "tmod_mat_.h"

// Subroutines in this file are machine-dependent and known to work on amd64
//  with FLINT ver 2.4.1 and 2.4.3

// Some code borrowed from flint nmod_mat/*.c, data layout is different

#include "tmod_mat_.h"

void 
tmod_mat_init(tmod_mat_t mat, long rows, long cols)
 {
  if(rows && cols)
   {
    slong i;
    mp_limb_t* e;
    mat->entries = flint_calloc(rows * cols, sizeof(mp_limb_t));
    mat->rows = flint_malloc(rows * sizeof(mp_limb_t *));
    e = mat->entries;
    for (i = 0; i < rows; i++, e += cols)
      mat->rows[i] = e;
   }
  else
   mat->entries = NULL;
  mat->r = rows;
  mat->c = cols;
 }

void 
tmod_mat_init_fast(tmod_mat_t mat, long rows, long cols)
// same result as tmod_mat_init(), only matrice is somewhat random rather than 
//  zero
 {
  slong i;
  mp_limb_t* e = mat->entries = flint_malloc(rows * cols * sizeof(mp_limb_t));
  mat->rows = flint_malloc( rows * sizeof(mp_limb_t*) );
  mat->r = rows;
  mat->c = cols;
  for(i=0; i<rows;  i++, e += cols)
    mat->rows[i] = e;
 }

void
tmod_mat_virginize(tmod_mat_t mat)
 {
  slong i,rows=mat->r,cols=mat->c;
  mp_limb_t* e = mat->entries;
  for(i=0; i<rows;  i++, e += cols)
    mat->rows[i] = e;
 }

void
tmod_mat_clear(tmod_mat_t mat)
 {
  if (mat->entries)
   {
    flint_free(mat->entries);
    flint_free(mat->rows);
   }
 }

mp_limb_t 
tmod_mat_diag_product_ZZ_ui(const tmod_mat_t m)
// counts product of diagonal entries
// m: count of rows >= count of columns > 0
 {
  slong cc=m->c-1;
  mp_limb_t** rows=m->rows;
  mp_limb_t r=rows[cc][cc];
  for( ; cc--; )
   r *= rows[cc][cc];
  return r;
 }

mp_limb_t
fmpz_to_t(const fmpz_t f)
/*
 This function calculates f modulo 2**64 on amd64. 
 Don't know what it does on other arch
*/
 {
  register slong n=(slong)(*f);
  register mp_limb_t m;
  if(!COEFF_IS_MPZ(n))
   {
    m=n;
    return m;
   }
  m=flint_mpz_get_ui(COEFF_TO_PTR(n));
  m *= (mp_limb_t)fmpz_sgn(f);
  return m;
 }

#include "tmod_mat_PLU.c"
#include "tmod_mat_window_unsh.c"
#include "tmod_mat_invert_LU.c"
