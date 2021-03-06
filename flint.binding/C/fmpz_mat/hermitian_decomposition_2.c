// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <string.h>
#if WANT_ASSERT_IN_HERM_DEC_2
 #include <assert.h>
 #define ASSERT(x) assert(x)
#else
 #define ASSERT(x)
#endif
#include <flint/flint.h>
#include "../nmod_mat/nmod_mat_.h"
#include "../fmpz_mat/fmpz_mat_.h"

#define LOUD_smallDet_matrice 0
#define LOUD_HERMITIAN_DECOMPOSITION 0

static __inline__ void
inverse_smallDet_HNF(fmpz_mat_t r,fmpz_t d,const nmod_mat_t a)
// r virgin
 {
  #if LOUD_smallDet_matrice
   flint_printf("small-det HNF dirty:\n");
   nmod_mat_print_pretty(a);
  #endif
  slong n=a->r,i,j;
  // copy diagonal and over-diagonal entries of a to r
  fmpz* rP=r->entries;
  mp_limb_t const* aP;
  for(i=0;i<n;i++)
   {
    // zero-out i entries, shift rP accordingly
    for(j=i;j--;rP++)
     fmpz_set_ui(rP, UWORD(0));
    aP=a->rows[i]+i;
    // copy n-1-i entries
    for(j=i;j<n;j++,aP++,rP++)
     fmpz_set_ui(rP, aP[0]);
   }
  #if LOUD_smallDet_matrice
   flint_printf("small-det HNF:\n");
   fmpz_mat_print_pretty(r);
  #endif
  fmpz_triU_inverse_smallDet(r,d, r);
  #if LOUD_smallDet_matrice
   flint_printf("it's inverse:\n");
   fmpz_mat_print_pretty(r);
  #endif
 }

static __inline__ void
_divide_away(fmpz_mat_t R,const fmpz_mat_t A,const fmpz_mat_t B,const fmpz_t D,
  slong n)
/*
 R, B virgin, n=dimension
 R = A*B/D, D fits 1 limb
*/
 {
  //fmpz_mat_mul() appears to be somewhat optimized, use it
  #if LOUD_smallDet_matrice
   flint_printf("multiplying:\n");
   fmpz_mat_print_pretty(A);
   flint_printf("by\n");
   fmpz_mat_print_pretty(B);
  #endif
  fmpz_mat_mul(R,A,B);
  #if LOUD_smallDet_matrice
   flint_printf("\nmultiply result:\n");
   fmpz_mat_print_pretty(R);
  #endif
  assert( fmpz_size(D)==1 );
  mp_limb_t d=fmpz_to_t(D );
  fmpz_mat_scalar_divexact_ui_2arg(R,d);
  #if LOUD_smallDet_matrice
   flint_printf("\ndivision result:\n");
   fmpz_mat_print_pretty(R);
  #endif
 }

int
fmpz_mat_hermitian_decomposition_2(fmpz_mat_t b,fmpz_t r, const fmpz_mat_t m)
/*
On entry b uninitialised.

Find hermitian decomposition of m with respect to 2**126: m=b*h.

If decomposition is trivial, 
 set b to shallow copy of m
 set r = 1
 return 0.

If det m is suspected to be divided by 2**127, 
 set b->r=0
 set r to a known divisor of det m (which happens to be 2**126)
 return 0

Otherwise 
 allocate b, 
 set b to m/h, 
 set r=det h
 return 1
*/
 {
  nmod_mat_t a; 
  nmod_mat_mod_t_half(a,m);
  #if LOUD_HERMITIAN_DECOMPOSITION
   gmp_printf("fmpz_mat_hermitian_decomposition_2():"
    " source matrice modulo T\n");
   nmod_mat_print_pretty(a);
  #endif
  // a must be cleared
  (void)nmod_mat_HNF(a);
  mp_limb_t d0=nmod_mat_diag_product_ZZ_ui(a);
  if(d0==1)
   {
    #if LOUD_HERMITIAN_DECOMPOSITION
     flint_printf(
      "fmpz_mat_hermitian_decomposition_2(): trivial case --- odd det\n");
    #endif
    memcpy(b, m, sizeof(fmpz_mat_struct));
    nmod_mat_clear(a);
    fmpz_set_ui(r,UWORD(1));
    return 0;
   }
  #if LOUD_HERMITIAN_DECOMPOSITION
   gmp_printf("fmpz_mat_hermitian_decomposition_2(): d0=%Mu factor a=\n",d0);
   nmod_mat_print_pretty(a);
  #endif
  slong n=a->r;
  fmpz_mat_t i; fmpz_mat_init(i,n,n);
  fmpz_t di; fmpz_init(di);
  inverse_smallDet_HNF(i,di, a);
  #if LOUD_HERMITIAN_DECOMPOSITION
   gmp_printf("fmpz_mat_hermitian_decomposition_2(): inv a=\n");
   nmod_mat_print_pretty(i);
  #endif
  nmod_mat_clear(a);
  // di, i must be cleared
  fmpz_mat_init(b,n,n);
  // b, di, i must be cleared
  _divide_away(b,m,i,di,n);
  if( d0 < UWORD(1)<<(FLINT_BITS-1) )
   {
    fmpz_clear(di);
    fmpz_mat_clear(i);
    fmpz_set_ui(r,d0);
    return 1;
   }
  nmod_mat_mod_t_half(a,b);
  // a, b, di, i must be cleared
  (void)nmod_mat_HNF(a);
  mp_limb_t d1=nmod_mat_diag_product_ZZ_ui(a);
  // r=d0*d1
  fmpz_set_ui(r,d0); fmpz_mul_ui(r,r,d1);
  if( d1 < UWORD(1)<<(FLINT_BITS-1) )
   {
    // b = b/a if d1>1
    if(d1>1)
     {
      inverse_smallDet_HNF(i,di, a);
      nmod_mat_clear(a);
      fmpz_mat_t c; fmpz_mat_init(c,n,n);
      _divide_away(c,b,i,di,n);
      fmpz_mat_clear(b);
      memcpy(b, c, sizeof(fmpz_mat_struct));
      fmpz_clear(di);
      fmpz_mat_clear(i);
     }
    else
     {
      nmod_mat_clear(a);
      fmpz_clear(di);
      fmpz_mat_clear(i);
     }
    return 1;
   }
  nmod_mat_clear(a);
  fmpz_mat_clear(b); b->r=0;
  fmpz_clear(di);
  fmpz_mat_clear(i);
  return 0;
 }

int nmod_mat_HNF_mod_t(nmod_mat_t A,fmpz_t d);

int
fmpz_mat_hermitian_decomposition_2_64(fmpz_mat_t b,fmpz_t r, const fmpz_mat_t m)
/*
On entry b uninitialised.

Attempt to decompose m into product b*h such that 
 det h is a degree of 2 
 det b is odd
 h entries fit 1 word

If diagonal entries of HNF of m modulo 2**64 are less than 2**64, decomposition 
 shall be found

If h is identity, 
 set b to shallow copy of m
 set r to 1
 return 0

If no decomposition found, 
 set b->r to 0
 set r to a known divisor of m
 return 0

Otherwise 
 allocate b
 set b to m/h
 set r=det h
 return 1
*/
 {
  nmod_mat_t a; 
  nmod_mat_mod_t(a,m);
  // a must be cleared
  int e=nmod_mat_HNF_mod_t(a,r);
  if(e)
   {
    b->r=0;
    nmod_mat_clear(a);
    return 0;
   }
  if(0==fmpz_cmp_ui(r, UWORD(1)))
   {
    memcpy(b, m, sizeof(fmpz_mat_struct));
    nmod_mat_clear(a);
    return 0;
   }
  slong n=a->r;
  fmpz_mat_t i; fmpz_mat_init(i,n,n);
  fmpz_t di; fmpz_init(di);
  inverse_smallDet_HNF(i,di, a);
  nmod_mat_clear(a);
  // di, i must be cleared
  fmpz_mat_init(b,n,n);
  _divide_away(b,m,i,di,n);
  fmpz_clear(di);
  fmpz_mat_clear(i);
  return 1;
 }

#undef NDEBUG
#undef ASSERT
