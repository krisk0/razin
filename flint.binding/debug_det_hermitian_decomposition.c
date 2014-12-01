// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// This program tests function fmpz_mat_det_hermitian_decomposition()
// Designed to run under valgrind

// To compile: gcc -O2 -march=native -lflint -lgmp flint_sage.so -g

#include <string.h>
#include "C/fmpz_mat/fmpz_mat_.h"
#define E_MAX 99

#define CALL_HD 1
#define CALL_FLINT 1

void fmpz_mat_det_hermitian_decomposition(fmpz_t r, const fmpz_mat_t a);

void
random_matr(fmpz_mat_t a,gmp_randstate_t rst,slong dim)
 {
  slong i,j;
  for(i=0;i<dim;i++)
   for(j=0;j<dim;j++)
    fmpz_set_si(fmpz_mat_entry(a,i,j), 
     (slong)(gmp_urandomm_ui(rst,2*E_MAX))-E_MAX);
 }

void
test_matr(fmpz_mat_t a)
 {
  fmpz_t z0; fmpz_init(z0);
  fmpz_t z1; fmpz_init(z1);
  #if CALL_FLINT
   fmpz_mat_det(z0,a);
  #endif
  #if CALL_HD
   fmpz_mat_det_hermitian_decomposition(z1,a);
  #endif
  #if CALL_HD+CALL_FLINT==2
   if( fmpz_cmp(z0,z1) )
    {
     flint_printf("det mismatch\n");
     abort();
    }
  #endif
  fmpz_clear(z1);
  fmpz_clear(z0);
 }

void 
do_with_dim(slong dim)
 {
  gmp_randstate_t rst; gmp_randinit_default(rst);
  gmp_randseed_ui(rst, 20141130-dim);
  fprintf(stderr,"dim=%d\n",dim);
  fmpz_mat_t a;
  fmpz_mat_init(a,dim,dim);
  slong i;
  for(i=4;i--;)
   {
    random_matr(a,rst,dim);
    test_matr(a);
   }
  fmpz_mat_clear(a);
  gmp_randclear(rst);
 }

int main()
 {
  slong i;
  for(i=3;i<501;i++)
   do_with_dim(i);
  flint_printf("test passed\n");
 }
