// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// This program tests procedure mpz_mod_2x()

// To compile: gcc -O2 -march=native -lflint -lgmp

#include "C/fmpz/fmpz_.h"

gmp_randstate_t rst;

void
do_test(mpz_t x,slong mu,mpz_t M)
 {
  mpz_t good; mpz_init(good);
  mpz_t baad; mpz_init(baad);
  mpz_mod(good,x,M);
  mpz_set(baad,x); mpz_mod_2x(baad,mu);
  if(mpz_cmp(good,baad))
   {
    gmp_printf("M=%ZX\n",M);
    gmp_printf("mu=%X x=%ZX good=%ZX baad=%ZX\n",mu,x,good,baad);
    abort();
   }
  mpz_clear(baad);
  mpz_clear(good);
 }

void
test(slong mu)
 {
  mpz_t M; mpz_init_set_ui(M,1); mpz_mul_2exp(M,M,mu);
  mpz_t a; mpz_init(a);
  
  do_test(a,mu,M);
  mpz_set_ui(a,1);
  do_test(a,mu,M);
  mpz_mul_2exp(a,a,mu-1);
  do_test(a,mu,M);
  mpz_sub_ui(a,M,1);  
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
  mpz_clear(M);
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
