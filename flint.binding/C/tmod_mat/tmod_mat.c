// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <gmp.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

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
// I am not a maniac, I am optimizer
 {
  slong i;
  mp_limb_t* e = mat->entries = flint_malloc(rows * cols * sizeof(mp_limb_t));
  mat->rows = flint_malloc( rows * sizeof(mp_limb_t*) );
  mat->r = rows;
  mat->c = cols;
  for (i = 0; i < rows; i++, e += cols)
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
This function calculates f modulo 2**64 on amd64, at least for FLINT ver 
 2.4.1 and 2.4.3. Don't know what it does on other arch
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

static mp_limb_t inv_mod__2_64__tab[13]={
 1,0xAAAAAAAAAAAAAAAB,0xCCCCCCCCCCCCCCCD,0x6DB6DB6DB6DB6DB7,0x8E38E38E38E38E39
 ,0x2E8BA2E8BA2E8BA3,0x4EC4EC4EC4EC4EC5,0xEEEEEEEEEEEEEEEF,0xF0F0F0F0F0F0F0F1,
 0x86BCA1AF286BCA1B,0xCF3CF3CF3CF3CF3D,0xD37A6F4DE9BD37A7,0x8F5C28F5C28F5C29};

static __inline__ mp_limb_t 
t_invmod(mp_limb_t a)
 {
  #define m 64
  mp_limb_t s;
  mp_limb_t t=a-1;
  if(t<26)
   return inv_mod__2_64__tab[ t>>1 ];
  // inspired by public domain code created by Jean-Guillaume Dumas
  count_trailing_zeros(s,t);
  t >>= s; 
  //assert( a == (1<<s)*t + 1 );
  //assert( t&1 );
  mp_limb_t U=2-a;
  mp_limb_t amone=a-1;
  long i,i_max;
  i_max=m/s;
  for(i=1;i<i_max;i <<= 1)
   {
    amone *= amone;
    U *= amone+1;
   }
  #undef m
  return U;
 }

#include "tmod_mat_PLU.c"
#include "tmod_mat_window_unsh.c"
#include "tmod_mat_invert_LU.c"
