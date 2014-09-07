// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef FMPZ__H
#define FMPZ__H

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

#endif
