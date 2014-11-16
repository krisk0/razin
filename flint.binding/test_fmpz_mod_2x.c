// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// This program tests procedure fmpz_mod_2x()

// To compile: gcc -O2 -march=native -lflint -lgmp

#include "C/fmpz/fmpz_.h"

gmp_randstate_t rst;

void
do_test(mpz_t y,slong mu,fmpz_t M)
 {
  fmpz_t good; fmpz_init(good);
  fmpz_t baad; fmpz_init(baad);
  fmpz_t x; fmpz_set_mpz(x,y);
  fmpz_mod(good,x,M);
  fmpz_set(baad,x); fmpz_mod_2x(baad,mu);
  if(fmpz_cmp(good,baad))
   {
    flint_printf("M="); fmpz_print(M); flint_printf("\n");
    gmp_printf("mu=%X x=%ZX good=%ZX baad=%ZX\n",mu,y,
     COEFF_TO_PTR(*good),
     COEFF_TO_PTR(*baad)
     );
    abort();
   }
  fmpz_clear(x);
  fmpz_clear(baad);
  fmpz_clear(good);
 }

void
test(slong mu)
 {
  fmpz_t M; fmpz_init_set_ui(M,1); fmpz_mul_2exp(M,M,mu);
  mpz_t a; mpz_init(a);
  
  do_test(a,mu,M);
  mpz_set_ui(a,1);
  do_test(a,mu,M);
  mpz_mul_2exp(a,a,mu-1);
  do_test(a,mu,M);
  fmpz_get_mpz(a,M);
  mpz_sub_ui(a,a,1);  
  do_test(a,mu,M);
  mpz_add_ui(a,a,1);  
  do_test(a,mu,M);
  mpz_add_ui(a,a,1);  
  do_test(a,mu,M);
  mpz_mul_ui(a,a,2);
  do_test(a,mu,M);
  if(mu>3)
   {
    slong i;
    for(i=10;i--;)
     {
      mpz_urandomb(a,rst,mu+1);
      do_test(a,mu,M);
     }
   }
  mpz_clear(a);
  fmpz_clear(M);
 }

int main()
 {
  gmp_randinit_default(rst);
  slong i;
  for(i=2;i<201;i++)
   test(i);
  flint_printf("Test passed\n");
  return 0;
 }
