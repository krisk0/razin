// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// This program tests subroutine fmpz_get_mpfr() and friends

#include <assert.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include "./C/fmpz/get_mpfr.c"
#include "./C/fmpz/set_mpfr.c"

flint_rand_t st;
fmpz_t two__64;

#define BITL(a) fmpz_sizeinbase(a,2)

#define IF_MORE_THAN_2                                    \
 if( fmpz_cmp_ui(m,UWORD(2)) > 0 )                          \
  {                                                          \
   flint_printf("n size=%d difference=",BITL(n));             \
   fmpz_print(m);                                              \
   flint_printf(" (%d)\n",BITL(m));                             \
  }

#define IF_MORE_THAN_2_AGAIN                               \
 if( fmpz_cmp_ui(m,UWORD(2)) > 0 )                           \
  {                                                            \
   flint_printf("|n| size=%d difference=",fmpz_sizeinbase(n,2)); \
   fmpz_print(m);                                               \
   flint_printf(" (%d)\n",BITL(m));                            \
  }

void 
test_1(fmpz_t n)
 {
  mpfr_t q0; mpfr_init(q0); 
  mpfr_t q1; mpfr_init(q1);
  fmpz_t m; fmpz_init(m);
  int c;
  
  fmpz_get_mpfr(     q0, n   , MPFR_RNDA );
  fmpz_get_mpfr_3arg(q1, n[0], MPFR_RNDA );
  assert( mpfr_equal_p( q0, q1 ) );
  fmpz_set_mpfr(m, q0, MPFR_RNDA);
  if( fmpz_cmp_ui( n, WORD(0) ) >= 0 )
   {
    if( fmpz_cmp(n,m) > 0 )
     {
      flint_printf("RNDA test failed, n="); fmpz_print(n);
      flint_printf(", m="); fmpz_print(m); flint_printf("\n");
      assert(0);
     }
    fmpz_sub(m, m, n);
    IF_MORE_THAN_2;
   }
  else 
   {
    if( fmpz_cmp(n,m) < 0 )
     {
      flint_printf("RNDA test failed, n="); fmpz_print(n);
      flint_printf(", m="); fmpz_print(m); flint_printf("\n");
      assert(0);
     }
    fmpz_sub(m, n, m);
    IF_MORE_THAN_2;
   }
   
  fmpz_get_mpfr( q0, n, MPFR_RNDZ );
  fmpz_set_mpfr( m, q0, MPFR_RNDZ );
  c=fmpz_cmp_si( n, WORD(0) );
  if(c==0)
   assert(0 == fmpz_cmp_si(m,WORD(0)));
  else
   {
    if(c<0)
     {
      fmpz_neg(n, n);
      fmpz_neg(m, m);
     }
    assert( fmpz_cmp(n,m) >= 0 );
    if(c)
     assert( fmpz_cmp_si(m,WORD(0)) >= 0 );
    fmpz_sub(m, n, m);
    IF_MORE_THAN_2_AGAIN;
   }
  
  fmpz_clear(m);
  mpfr_clear(q1);
  mpfr_clear(q0);
 }

void 
grow(fmpz_t x)
 {
  fmpz_mul(x, x, two__64);
  fmpz_add_ui(x, x, n_randlimb(st) );
 }

void 
test_0( slong x )
 {
  fmpz_t t; fmpz_init(t); fmpz_set_si( t, x );
  test_1( t );
  grow( t );
  test_1( t );
  grow( t );
  test_1( t );
  fmpz_clear(t);
 }

slong 
R()
 {
  mp_limb_t t;
  slong r;
  while(1)
   {
    t=n_randlimb(st);
    if(t<COEFF_MAX)
     break;
   }
  r=t;
  if( n_randlimb(st) & 1 )
   r=-r;
  return r;
 }

int 
main()
 {
  flint_randinit(st);
  fmpz_init_set_ui( two__64, UWORD(1)<<32 ); 
  fmpz_mul( two__64, two__64, two__64 );
  int i;
  test_0( 0 );
  for(i=100;i--;)
   test_0( R() );
  fmpz_clear(two__64);
  flint_printf("Test passed\n");
  return 0;
 }
