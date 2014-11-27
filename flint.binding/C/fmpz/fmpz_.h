// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef FMPZ__H
#define FMPZ__H

// mpz/mpq/fmpz/mpfr manipulation macros/functions

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
decrease_bound_mpz(mpfr_t b,mpfr_prec_t pr,mpz_t d)
 {
  mpfr_t dF; mpfr_init2(dF, mpz_sizeinbase(d,2)+3 );
  mpfr_t log2_d; mpfr_init2(log2_d,pr);

  mpfr_set_z(dF, d, MPFR_RNDZ);
  mpfr_log2(log2_d, dF, MPFR_RNDZ);
  mpfr_sub(b, b, log2_d, MPFR_RNDU);

  mpfr_clear(log2_d);
  mpfr_clear(dF);
 }

static __inline__ void
decrease_bound_mpz_2arg(mpfr_t b,mpz_t d)
 {
  mpfr_t dF; mpfr_init2(dF, mpz_sizeinbase(d,2)+3 );
  mpfr_t log2_d; mpfr_init2(log2_d,mpfr_get_prec(b));

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

#define SURELY_REDUCE_SIZE \
 while( r_limbs-- && (0==r->_mp_d[r_limbs]) ) \
  ;                                            \
 if(r_limbs >= 0)                               \
  r->_mp_size=r_limbs+1;                        \
 else                                           \
  r->_mp_size=0;                                

__inline__ void
mpz_mod_2x(mpz_t r,slong x)
// reduce positive r modulo 2**x
 {
  slong r_limbs=r->_mp_size;
  slong x_limbs=(x+FLINT_BITS-1)/FLINT_BITS;
  if(r_limbs<x_limbs)
   return;
  if(r_limbs>x_limbs)
   r_limbs=r->_mp_size=x_limbs;
  // reduce senior limb
  x_limbs=FLINT_BITS*x_limbs-x; // this many senior bits must be cut-off
  if(x_limbs)
   {
    mp_limb_t m=( UWORD(1)<<( FLINT_BITS - x_limbs ))-1;
    mp_limb_t w=r->_mp_d[--r_limbs];
    mp_limb_t u=w&m;
    if( u<w )
     {
      r->_mp_d[r_limbs]=u;
      if(u)
       return;
      // must decrease r->_mp_size that now equals r_limbs+1
      SURELY_REDUCE_SIZE
      return;
     }
    if(u)
     return;
    SURELY_REDUCE_SIZE
   }
  else
   {
    // maybe decrease r->_mp_size that now equals r_limbs
    SURELY_REDUCE_SIZE
   }
 }

#undef SURELY_REDUCE_SIZE

__inline__ void
fmpz_mod_2x(fmpz_t r,slong x)
// reduce positive r modulo 2**x
 {
  fmpz c=*r;
  if( COEFF_IS_MPZ(c) )
   {
    mpz_mod_2x(COEFF_TO_PTR(c), x);
    // r might get very small
    _fmpz_demote_val(r);
   }
  else
   {
    if(x >= FLINT_BITS)
     return;   // r already smaller than 2**64<=2**x
    *r &= (UWORD(1)<<x)-1;
   }
 }

#define MPFR_INIT(x) { mpfr_init(x); mpfr_set_ui(x,0,MPFR_RNDU); }
#define MPFR_INIT2(x,y) { mpfr_init2(x,y); mpfr_set_ui(x,0,MPFR_RNDU); }

__inline__ mp_limb_t
mpz_to_t(mpz_ptr x)
 {
  slong s=x->_mp_size;
  if( 0==s )
   return 0;
  return s>0 ? x->_mp_d[0] : -x->_mp_d[0];
 }

__inline__ int
is_identity_permutation(mp_limb_t* p,slong n)
 {
  slong i;
  for(i=0; i<n; i++)
   {
    if( p[i] != (mp_limb_t)i )
     return 0;
   }
  return 1;
 }

__inline__ int
count_permutation_parity(mp_limb_t* P,slong n)
/*
linear in space and time

inspired by public-domain code found at 
 http://www.cap-lore.com/code/ocaml/parity.html

P represented with lower-line
return 1 iff P is even, -1 otherwise
*/
 {
  int p=0; 
  mp_limb_t* v=flint_calloc(n,sizeof(mp_limb_t));
  slong j=n; 
  while(j--)
   if(v[j])
    ++p;
   else
    {
     slong x=j;
     do 
      {
       x = P[x]; v[x]=1;
      }
     while(x!=j);
    }
  flint_free(v);
  return 1-2*(p&1);
 }

__inline__ void
muladd_mpz_fmpz(mpz_t r,mpz_t a,fmpz_t b)
// r += a*b
 {
  fmpz bb=*b;
  if( COEFF_IS_MPZ(bb) )
   mpz_addmul(r, a, COEFF_TO_PTR(bb));
  else
   {
    if(bb<0)
     flint_mpz_submul_ui(r,a,-bb);
    else
     flint_mpz_addmul_ui(r,a,bb);
   }
 }

mp_limb_t
fmpz_to_t(const fmpz_t f)
/*
 This function calculates f modulo 2**FLINT_BITS on amd64. 
 Don't know what it does on other arch
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

void __inline__ static
MulMod_2x_stupid(mpz_t tgt,const mpz_t sou,slong log2_M)
// multiply modulo 2**log2_M
 {
  mpz_mul(tgt,tgt,sou);
  mpz_mod_2x(tgt,log2_M);
 }

/*
MulMod_2x_positive(tgt,sou,log2_M):
 tgt is zero or positive and has enough place for log2_M bits
 sou is positive and has enough place for log2_M bits

auto-generated file defining subroutine MulMod_2x_positive():
*/
#include "MulMod_2x_positive.c"

#endif
