// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/
// Copyright Денис Крыськов 2014

// This file is based on det_modular_given_divisor.c which is part of FLINT and
//  owned by F.Johansson

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include "../ulong_extras/ulong_extras_.h"
#include <assert.h>
void nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim);
mp_limb_t nmod_mat_det_mod_pk(nmod_mat_t M,const p_k_pk_t pp,mp_limb_t* scrtch);

static __inline__ void
max_degree(p_k_pk_t s)
 {
  mp_limb_t b=UWORD_MAX / s.p;
  s.k=1;
  s.p_deg_k=s.p;
  while( s.p_deg_k <= b )
   {
    s.p_deg_k *= s.p;
    ++s.k;
   }
  assert( 0 == s.p_deg_k % s.p );
 }

static __inline__ mp_limb_t
select_prime_and_degree(p_k_pk_t pp,nmod_t mod,const fmpz_t divisor)
 {
  mp_limb_t r,t;
  while(1)
   {
    pp.p = n_nextprime(pp.p, 0); // in FLINT 2.4.4 this is prime
    max_degree( pp );
    r=fmpz_fdiv_ui( divisor, pp.p_deg_k );
    if(r % pp.p)
     {
      count_leading_zeros( t, pp.p_deg_k );
      mod.n = pp.p_deg_k << t;
      mod.ninv = n_preinvert_limb(mod.n);
      return inv_mod_pk_3arg(r,pp,mod);
     }
   }
 }

void 
fmpz_mat_det_modular_given_divisor_3arg(fmpz_t det,const fmpz_mat_t A,
  fmpz_t divisor)
 {
  fmpz_t bound, prod, x, xnew;
  mp_limb_t xmod;
  mp_limb_t divisor_inv; // stands for n_invmod(fmpz_fdiv_ui(d, p)
  p_k_pk_t pp; pp.p=2;
  nmod_mat_t Amod; Amod->mod.norm=0;
  slong dim = A->r;
  mp_limb_t* scratch=flint_malloc( 4*(dim-4)*sizeof(mp_limb_t) );
  fmpz_init(bound);
  // TODO: check if line below brings in useless branch
  fmpz_init_set_ui( prod, UWORD(1) );
  fmpz_init(x);
  fmpz_init(xnew);
  nmod_mat_init_square_2arg(Amod,dim);
  
  /* bound = 2 * max abs det(A) / divisor */
  // TODO: use sum of logarithms instead of big number product
  fmpz_mat_det_bound(bound, A);
  fmpz_mul_ui(bound, bound, UWORD(2));
  fmpz_cdiv_q(bound, bound, divisor);

  while(fmpz_cmp(prod, bound) <= 0)
   {
    divisor_inv=select_prime_and_degree( pp, Amod->mod, divisor );
    fmpz_mat_get_nmod_mat(Amod, A);
    xmod=nmod_mat_det_mod_pk(Amod,pp,scratch);
    xmod=n_mulmod2_preinv(xmod,divisor_inv, Amod->mod.n,Amod->mod.ninv);
    
    fmpz_CRT_ui(xnew, x, prod, xmod, pp.p_deg_k, 1);
    fmpz_mul_ui(prod, prod, pp.p_deg_k);
    fmpz_set(x, xnew);
   }
  
  nmod_mat_clear(Amod);
  fmpz_clear(xnew);
  fmpz_clear(x);
  fmpz_clear(prod);
  fmpz_clear(bound);
  flint_free(scratch);
 }
