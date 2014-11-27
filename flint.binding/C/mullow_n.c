// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <stdlib.h>
#include <flint/flint.h>

/*
 attempt to multiply two numbers modulo 2**128 with undocumented GMP 
  subroutine mpn_mullow_n... no, mpn_mullo_n ... no, __gmpn_mullo_n

 return 0 if the trick worked correctly
*/

#define FUNC __gmpn_mullo_n
void FUNC(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n);
#define PTR(x) ((x)->_mp_d)

int main(int ignored,char** v)
 {
  mpz_t a,b,c;
  mpz_init2(a,2*FLINT_BITS); mpz_set_ui(a, atol(v[1]));
  mpz_init(b);
  mpz_init_set_ui(c,1); mpz_mul_2exp(c,c,2*FLINT_BITS);
  mpz_invert(b,a,c); // a*b must be 1 modulo c
  mpz_sub_ui(c,c,1);
  mp_limb_t* c_inside=PTR(c);
  FUNC(c_inside, PTR(b), PTR(a), 2);
  if( 0==c_inside[1] )
   {
    --c->_mp_size;
    if( 0==c_inside[0] )
     --c->_mp_size;
   }
  return mpz_cmp_ui(c,1);
  // memory leak: not calling mpz_clear()
 }
