// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include "tmod_mat_.h"
#include "../fmpz/fmpz_.h"

// Subroutines in this file are machine-dependent and known to work on amd64

// Some code borrowed from flint nmod_mat/*.c, data layout is different


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
tmod_mat_init_set_fmpz_mat(tmod_mat_t tgt, const fmpz_mat sou)
 {
  slong i, j, rows=sou->r, cols=sou->c;
  mp_limb_t* e = tgt->entries = flint_malloc(rows * cols * sizeof(mp_limb_t));
  mp_limb_t* tgt_rows = tgt->rows = flint_malloc( rows * sizeof(mp_limb_t*) );
  tgt->r = rows;
  tgt->c = cols;
  fmpz* souP;
  for(i=0; i<rows;  i++, e += cols)
   {
    tgt_rows[i] = e;
    for(j=0,souP=sou->rows[i]; j<cols; j++,souP++)
     e[j]=fmpz_to_t( souP );
   }
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
dot_modulo_t_kind0( mp_limb_t*a, fmpz* b, slong n )
// a * b
 {
  mp_limb_t s=0;
  for(;n--; a++,b++)
   s += a[0] * fmpz_to_t( b );
  return s;
 }

#include "tmod_mat_PLU.c"
#include "tmod_mat_window_unsh.c"
#include "tmod_mat_invert_LU.c"
