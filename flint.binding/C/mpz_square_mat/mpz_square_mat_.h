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
  slong* mark;
 }
mpz_square_mat_struct;

typedef mpz_square_mat_struct mpz_square_mat_t[1];

static __inline__ void 
mpz_square_mat_lazy_init_set(mpz_square_mat_t R,const fmpz_mat_t S)
 {
  R->mark=0;
  slong i,j,siz=R->r=S->r;
  mpz_ptr e=flint_malloc(siz*siz*sizeof(__mpz_struct));
  R->entries=e;
  R->rows=flint_malloc(siz*sizeof(mpz_ptr));
  fmpz* g;
  for(i=0;i<siz;i++)
   {
    R->rows[i]=e;
    g=S->rows[i];
    for(j=siz;j--;e++,g++)
     {
      mpz_init( e );
      fmpz_get_mpz( e, g );
     }
   }
 }

static __inline__ void 
mpz_square_mat_clear(mpz_square_mat_t m)
 {
  slong i=m->r; ;
  for(i=i*i; i--; )
   mpz_clear(m->entries+i);
  flint_free(m->entries);
  flint_free(m->rows);
  flint_free(m->mark);
 }

static __inline__ void 
mpz_square_mat_init_transpose(mpz_square_mat_t R,const mpz_square_mat_t S)
 {
  R->mark=0;
  slong i,j,siz=R->r=S->r;
  mpz_ptr e=flint_malloc(siz*siz*sizeof(__mpz_struct));
  R->entries=e;
  R->rows=flint_malloc(siz*sizeof(mpz_ptr));
  mpz_ptr* const s_rows=S->rows;
  mpz_ptr g=s_rows[0];
  // transfer column 0, fill R->rows[]
  for(i=0;i<siz;i++,e += siz,g++)
   {
    R->rows[i]=e;
    mpz_init_set(e, g);
   }
  for(j=1;j<siz;j++)
   {
    // transfer column j
    g=s_rows[j];
    e=R->entries+j;
    for(i=siz;i--;e += siz,g++)
     mpz_init_set(e, g);
   }
 }

static __inline__ void
mpz_square_mat_mark_biggest(mpz_square_mat_t m)
// mark[i] := FLINT_BITS*(2+limb-size of biggest entry in row i)
 {
  slong dim=m->r,i,j;
  m->mark=flint_malloc( sizeof(slong)*dim );
  mpz_ptr p;
  slong b=0,c;
  for(i=0;i<dim;i++)
   {
    for(j=dim,p=m->rows[i];j--;)
     {
      c=mpz_size( p++ );
      if(c>b)
       b=c;
     }
    m->mark[i]=FLINT_BITS*(2+b); // TODO: maybe add 3 instead of 2?
   }
 }

static __inline__ void
mpz_square_mat_negate(mpz_square_mat_t m)
 {
  slong s=m->r; s=s*s;
  mpz_ptr e=m->entries;
  for(;s--;e++)
   mpz_neg( e, e );
 }

void mpz_square_mat_mul_vec_mat_modulo(mpz_ptr t,
  const mpz_ptr v,const mpz_square_mat_t A,const mpz_t m);

// macro used by Dixon lifting subroutines
__inline__ static void
det_divisor_init_b(mpz_ptr b,slong n,mpz_square_mat_t m)
 {
  slong i;
  mpz_ptr t;
  slong* s=m->mark;
  for(i=0,t=b;i<n;i++,t++)
   {
    mpz_init2(t, s[i]);
    mpz_set_si(t,  2*(i&1)-1);
   }
 }

#endif
