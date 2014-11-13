// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef FMPZ__H
#define FMPZ__H

// mpz/mpq/fmpz manipulation macros/functions

#include <flint/flint.h>
#include <flint/fmpz.h>

// s should be fmpz, s evaluated multiple times
#define fmpz_get_mpfr_slave(r, s, rnd)    \
 if( COEFF_IS_MPZ(s) )                      \
   (void)mpfr_set_z(r, COEFF_TO_PTR(s), rnd); \
 else                                         \
   mpfr_set_si(r, s, rnd);

// _s should be fmpz, _s evaluated once
#define fmpz_get_mpfr_macro(r, _s, rnd) \
 {                                     \
  register fmpz _sR=_s;               \
  fmpz_get_mpfr_slave(r,_sR,rnd);    \
 }

void __inline__ static
mul_mpz_fmpz(mpz_t r,fmpz_t s)
// this trick is documented in FLINT manual
 {
  __mpz_struct *z;
  z=_fmpz_promote_val(s);
  mpz_mul(r,r,z);
  _fmpz_demote_val(s);
 }

void __inline__ static
fmpz_mod_2arg(fmpz_t tgt,const mpz_t m)
 {
  if(fmpz_fits_si(tgt))
   {
    mpz_t t; mpz_init_set_si(t,fmpz_get_si(tgt));
    mpz_mod(t,t,m);
    fmpz_set_mpz(tgt,t);
    mpz_clear(t);
   }
  else
   {
    __mpz_struct *t = _fmpz_promote(tgt);
    mpz_mod(t, t, m);
    _fmpz_demote_val(tgt);
   }
 }
 
int __inline__ static
abs_x_is_degree_of_2(const mpz_t x)
 {
  size_t s=mpz_size(x)-1;
  slong i;
  for(i=0;i<s;i++)
   if( mpz_getlimbn(x,i) )
    return 0;
  mp_limb_t y=mpz_getlimbn(x,s);
  // |x| is degree of 2 iff subtracting 1 from y decreases its bit-length
  return !( y&(y-1));
 }

void __inline__ static
fmpz_hex_print(const char* m,const fmpz_t n,int CR)
 {
  flint_printf(m);
  if(fmpz_fits_si(n))
   {
    // could use WORD_FMT like that: printf(WORD_FMT"X",small_num);
    mpz_t t; mpz_init_set_si(t,fmpz_get_si(n));
    gmp_printf("%ZX (%d)",t,mpz_sizeinbase(t,2));
    mpz_clear(t);
   }
  else
   {
    #define nn (fmpz*)n
    __mpz_struct* t = _fmpz_promote(nn);
    gmp_printf("%ZX (%d)",t,mpz_sizeinbase(t,2));
    _fmpz_demote_val(nn);
    #undef nn
   }
  if(CR)
   flint_printf("\n");
 }

void __inline__ static
mpz_hex_print(const char* m,const mpz_t n,int CR)
 {
  flint_printf(m);
  gmp_printf("%ZX (%d)",n,mpz_sizeinbase(n,2));
  if(CR)
   flint_printf("\n");
 }
 
ulong __inline__ static
fmpz_fdiv_ui_positive(const fmpz_t g, ulong h)
 {
  fmpz c1 = *g;
  if (!COEFF_IS_MPZ(c1))      // g is small 
   return c1 % h;
  else                        // g is large
   // should work correctly under Linux and Windows/MPIR
   return mpz_fdiv_ui(COEFF_TO_PTR(c1), h);
 }

int __inline__ static
cmp_positive_log2(const fmpz_t z,mp_limb_t y)
// both numbers positive
 {
  register mp_limb_t xx=(mp_limb_t)(*z);
  int lz;
  if(!COEFF_IS_MPZ(xx))
   {
    if(y>=FLINT_BITS)
     return -1;          // |z|<2**64, 2**y>2**64
    count_leading_zeros(lz,xx);
    lz=FLINT_BITS-lz;    // z bit-length = lz, hence z in 2**(lz-1)..2**lz-1
    if(lz <= y)
     return -1;
    if(lz-1 > y)
     return 1;
    // return 0 if xx is degree of 2, else return positive
    return ( xx & (xx-1) )>0;
   }
  mpz_ptr x=COEFF_TO_PTR( (slong)xx );
  xx=mpz_size(x);      // xx>1
  if(xx*FLINT_BITS <= y)
   return -1;          // z is too short, no need to further compare
  --xx;
  if( xx*FLINT_BITS > y )
   return 1;           // even if senior limb is 1, x is larger
  count_leading_zeros(lz,mpz_getlimbn(x,xx));
  xx=(xx*FLINT_BITS)+FLINT_BITS-lz;   // xx=|z| bit_length
  if(xx <= y)
   return -1;     // |z| fits y bits
  if(xx-1 > y)
   return 1;      // |z| > 2**y because its bit-length is y+2 or more
  // Now either |z| is degree of 2, and result should be 0
  //  or |z| is not degree of 2, and result should be 1
  return !abs_x_is_degree_of_2(x);
 }

static __inline__ void
decrease_bound_fmpz(mpfr_t b,mpfr_prec_t pr,mpz_t d)
 {
  mpfr_t dF; mpfr_init2(dF, mpz_sizeinbase(d,2));
  mpfr_t log2_d; mpfr_init2(log2_d,pr);

  mpfr_set_z(dF, d, MPFR_RNDZ);
  mpfr_log2(log2_d, dF, MPFR_RNDZ);
  mpfr_sub(b, b, log2_d, MPFR_RNDU);

  mpfr_clear(log2_d);
  mpfr_clear(dF);
 }

// small macro used by fmpz_mat_det_modular_given_divisor_8arg() and
//  fmpz_mat_det_modular_given_divisor_4arg()
static __inline__ void
mpz_fmpz_mul_det_2arg(mpz_t z,const fmpz_t x)
 {
  register slong xx=(slong)(*x);
  if(!COEFF_IS_MPZ(xx))
   {                              // x is small
    if(xx != WORD(1))
     mpz_mul_si(z,z,xx);
   }
  else
   mpz_mul(z,z,COEFF_TO_PTR(xx)); // x is big
 }

// macro used multiple times by *det_divisor*()
__inline__ static void
clear_mpz_array(mpz_ptr b,slong n)
 {
  slong i;
  mpz_ptr t;
  for(i=n,t=b;i--;t++)
   mpz_clear(t);
  flint_free(b);
 }

__inline__ static void
clear_mpfr_array(__mpfr_struct* b,slong n)
 {
  slong i;
  __mpfr_struct* t;
  for(i=n,t=b;i--;t++)
   mpfr_clear(t);
  flint_free(b);
 }

static __inline__ void
mpfr_copy_bound(mpfr_t r,const mpfr_t s)
 {
  mpfr_init2(r, mpfr_get_prec(s));
  mpfr_set(r, s, MPFR_RNDU);
 }

__inline__ static void
square_L2_fmpz(fmpz_t r,const fmpz* vec,slong n)
// on entry r=0. on exit r=squared L2 norm of the vector 
 {
  fmpz_mul(r,vec,vec);
  const fmpz* u;
  fmpz_t x; fmpz_init(x);
  for(u=vec+1,--n;n--;u++)
   {
    fmpz_mul(x,u,u);
    fmpz_add(r,r,x);
   }
  fmpz_clear(x);
 }

__inline__ static void
mpz_shift_right_1limb(mpz_t r)
 {
  slong limbs=r->_mp_size;
  if(0==limbs)
   return;
  if(limbs>0)
   {
    --limbs;
    --r->_mp_size;
   }
  else
   {
    limbs=(-limbs)-1;
    ++r->_mp_size;
   }
  mp_limb_t* n=r->_mp_d;
  for(;limbs--;++n)
   n[0]=n[1];
 }

__inline__ static mp_limb_t
mpz_mod_T(mpz_t r)
 {
  if( 0==r->_mp_size )
   return 0;
  mp_limb_t m;
  if( 0<r->_mp_size )
   m=r->_mp_d[0];
  else
   m=-r->_mp_d[0];
  return m;
 }

#endif
