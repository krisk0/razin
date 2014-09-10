// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz_mat.h>
#include <mpfr.h>
#include "../ulong_extras/ulong_extras_.h"
#include "../fmpz/fmpz_.h"
#include "../nmod_mat/nmod_mat_.h"
#undef NDEBUG
#include <assert.h>

static __inline__ void
decrease_bound_ui(mpfr_t b,mpfr_prec_t pr,mp_limb_t d)
 {
  mpfr_t dF; mpfr_init2(dF, FLINT_BITS);
  mpfr_t log2_d; mpfr_init2(log2_d,pr);

  mpfr_set_uj(dF, d, MPFR_RNDZ);
  mpfr_log2(log2_d, dF, MPFR_RNDZ);
  mpfr_sub(b, b, log2_d, MPFR_RNDU);

  mpfr_clear(log2_d);
  mpfr_clear(dF);
 }

static __inline__ int
comp_bound_ui(mpfr_t b,mp_limb_t q)
// 1 iff log2(q) rounded down <= b
 {
  // TODO: use MPFR_DECL_INIT instead of mpfr_init2
  mpfr_t log2; mpfr_init(log2);
  mpfr_t qF; mpfr_init2(qF, FLINT_BITS);
  mpfr_set_uj(qF, q, MPFR_RNDZ);
  mpfr_log2(log2, qF, MPFR_RNDZ);
  int r=(mpfr_cmp(log2,b)<0);
  mpfr_clear(qF);
  mpfr_clear(log2);
  return r;
 }

static __inline__ mp_limb_t
invert_det_divisor_modulo_pk(mpz_t dd,p_k_pk_t const* pp,nmod_t const* mod)
// take dd modulo p_deg_k, then invert it
 {
  mp_limb_t r,m=pp->p_deg_k;
  r=mpz_fdiv_ui(dd, m);
  if(pp->k==1)
   return n_invmod(r,m);
  else
   return inv_mod_pk_3arg(r,pp[0],mod[0]);
 }

static __inline__ void
mpz_fmpz_mul_2arg(mpz_t z,const fmpz_t x)
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

static __inline__ mp_limb_t
choose_prime_and_degree(p_k_pk_t* pp,nmod_t* mod,n_primes_rev_t it,
  const mpz_t divisor)
 {
  mp_limb_t r,t,r_mod_p;
  for(;;)
   {
    if( (pp->p = n_primes_rev_next(it)) == 1 )
     {
      flint_printf("Exception (choose_prime_and_degree): "
                     "Prime set exhausted\n");
      abort();
     }
    if( pp->p >= (UWORD(1)<<(FLINT_BITS/2)) )
     {
      pp->k=1;
      r=r_mod_p=mpz_fdiv_ui( divisor, pp->p_deg_k=pp->p );
     }
    else
     {
      max_degree( pp );
      r=mpz_fdiv_ui( divisor, pp->p_deg_k );
      r_mod_p=r % pp->p;
     }
    if(r_mod_p)
     {
      count_leading_zeros( t, pp->p_deg_k );
      mod->n = pp->p_deg_k << t;
      invert_limb(mod->ninv, mod->n);
      #if SPEEDUP_NMOD_RED3
       t = - mod->n;
       mod->norm = n_mulmod_preinv_4arg( t,t, mod->n,mod->ninv );
      #endif
      return inv_mod_pk_4arg(r,r_mod_p,pp[0],mod[0]);
     }
   }
 }

void 
fmpz_mat_det_modular_given_divisor_8arg(mpz_t det, nmod_mat_t Amod,
  mpfr_t hadamard_log2, mpfr_prec_t pr, p_k_pk_t* pp, n_primes_rev_t iT,
  mp_limb_t xmod, const fmpz_mat_t A)
/*
act like fmpz_mat_det_modular_given_divisor_4block(), but
 * decrease primes using iT, rather than increase them
 * don't count H.B., use hadamard_log2 which is 1+log2(H.B.)
 * sum logarithms instead of multiplying together found primes
 * re-use found prime pp->p and xmod which is determinant of A modulo pp->p

hadamard_log2 on entry is upper bound on log2(2*H.B)
              on exit, decreased to an unspecified value 

iT on entry just found prime pp->p
   possibly gets shifted 
*/
 {
  // loop bound = 2*H.B / known det divisor
  decrease_bound_fmpz(hadamard_log2,pr,det);

  // re-use known det A modulo pp->p_deg_k
  mp_limb_t divisor_inv=invert_det_divisor_modulo_pk(det,pp,&Amod->mod);
  xmod=n_mulmod2_preinv(xmod,divisor_inv, Amod->mod.n,Amod->mod.ninv);
  fmpz_t xnew,x;
  fmpz_init(xnew); 
  fmpz_init(x); 
  fmpz_t prod; fmpz_init_set_ui( prod, UWORD(1) );

  fmpz_CRT_ui(xnew, x, prod, xmod, pp->p_deg_k, 1);
  fmpz_set_ui(prod, pp->p_deg_k);
  fmpz_set(x, xnew);

  #if LOUD_DET_BOUND
   mpfr_printf("fmpz_mat_det_modular_given_divisor_8arg(): log2 bound=%Rf\n",
    hadamard_log2);
   slong primes_used=1;
  #endif
  // for orthogonal matrice the bound might be reached at this point.
  //  Attempt to skip main loop
  if(comp_bound_ui(hadamard_log2,pp->p_deg_k))
   {
    mp_limb_t* scratch=flint_malloc( 4*(A->r-4)*sizeof(mp_limb_t) );
    mp_limb_t bound=mpfr_get_uj(hadamard_log2,MPFR_RNDU);

    for(;;)
     {
      divisor_inv=choose_prime_and_degree( pp, &Amod->mod, iT, det );
      // TODO: optimize fmpz_mat_get_nmod_mat()
      fmpz_mat_get_nmod_mat(Amod, A);
      // TODO: call a faster subroutine instead of nmod_mat_det_mod_pk_4block()
      //  when pp->p is 64 bit long
      xmod=nmod_mat_det_mod_pk_4block(Amod,pp[0],scratch);
      xmod=n_mulmod2_preinv(xmod,divisor_inv, Amod->mod.n,Amod->mod.ninv);
      // TODO: rewrite fmpz_CRT_ui() -> mpz_CRT_ui_5arg()
      fmpz_CRT_ui(xnew, x, prod, xmod, pp->p_deg_k, 1);
      fmpz_mul_ui(prod, prod, pp->p_deg_k);
      #if LOUD_DET_BOUND
       primes_used++;
      #endif
      if(cmp_positive_log2(prod,bound) >= 0)
       break;
      fmpz_set(x, xnew);
     }
  
    flint_free(scratch);
   }

  #if LOUD_DET_BOUND
   flint_printf("fmpz_mat_det_modular_given_divisor_8arg() primes used: %d\n\n\n",
    primes_used);
  #endif
  fmpz_clear(prod);
  mpz_fmpz_mul_2arg(det,xnew);
  fmpz_clear(prod);
  fmpz_clear(x);
  fmpz_clear(xnew);
 }

#undef NDEBUG
