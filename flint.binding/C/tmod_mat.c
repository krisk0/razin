// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

//#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include "flint/ulong_extras.h"
#include "flint/fmpz.h"

typedef struct
 {
  mp_limb_t* entries;
  slong r;
  slong c;
  mp_limb_t ** rows;
 }
tmod_mat_struct;
typedef tmod_mat_struct tmod_mat_t[1];
#define tmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])

// Some code borrowed from flint-2.4.1/nmod_mat/

void 
tmod_mat_init(tmod_mat_t mat, long rows, long cols)
 {
  if(rows && cols)
   {
    long i;
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
tmod_mat_clear(tmod_mat_t mat)
 {
  if (mat->entries)
   {
    flint_free(mat->entries);
    flint_free(mat->rows);
   }
 }

mp_limb_t
fmpz_to_t(const fmpz_t f)
// This function calculates f modulo 2**64 on amd64. Don't know what it does on
//  other arch
 {
  register long n=(long)(*f);
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

mp_limb_t
t_invmod(mp_limb_t x)
 // returns x inverted modulo 2**64. x must be less than 2**63 and odd
 {
  mp_limb_t r=n_invmod(x,0x8000000000000000);
  if(r*x == 1)
   return r;
  return r*0x8000000000000001;
 }

#include "tmod_mat_PLU.c"
