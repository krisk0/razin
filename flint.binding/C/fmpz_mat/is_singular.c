// This program is part of RAZIN
// Copyright Денис Крыськов 2014

/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

// This file contains modified code of det_divisor.c which is part of FLINT and
//  owned by F.Johansson

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

#if GMP_LIMB_BITS == 64
#include "../ulong_extras/ulong_extras_.h"

void nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim);

slong __inline__ static
_call_solve_Dixon(const fmpz_mat_t A)
 {
  slong i,n=A->r;
  int success;
  fmpz_mat_t B,X;
  fmpz_mat_init(B, n, 1); fmpz_set_ui( B->entries, UWORD(1) );
  fmpz_mat_init(X, n, 1);
  fmpz_t mod; fmpz_init(mod);
  
  success = fmpz_mat_solve_dixon(X, mod, A, B);
  
  fmpz_clear(mod);
  fmpz_mat_clear(X);
  fmpz_mat_clear(B);
  return !success;
 }

slong 
fmpz_mat_is_singular(const fmpz_mat_t A)
 {
  slong dim=A->r,i;
  int j=1;
  if(dim<5)
   {
    fmpz_t d; fmpz_init(d);
    fmpz_mat_det_cofactor(d, A);
    j=!fmpz_cmp_ui(d, UWORD(0));
    fmpz_clear(d);
    return j;
   }
  mp_limb_t* scratch=flint_malloc( 4*(dim-4)*sizeof(mp_limb_t) );
  nmod_mat_t Amod; 
  p_k_pk_t pp; pp.k=1;
  #if !SPEEDUP_NMOD_RED3
   Amod->mod.norm=0;
  #endif
  nmod_mat_init_square_2arg(Amod,dim);
  mp_limb_t pr[10]={0x95,0x9F,0xBD,0xE9,0xFF,0x143,0x14D,0x1A1,0x1AD,0x1C5};
  mp_limb_t t;
  // test 10 biggest primes
  for(i=10;i--;)
   {
    pp.p=pp.p_deg_k=Amod->mod.n=pr[i] ^ UWORD(0xFFFFFFFFFFFFFE00);
    invert_limb(Amod->mod.ninv, pp.p);
    #if SPEEDUP_NMOD_RED3
     t = - pp.p;
     Amod->mod.norm = n_mulmod_preinv_4arg( t,t, pp.p,Amod->mod.ninv );
    #endif
    fmpz_mat_get_nmod_mat(Amod, A);
    if( nmod_mat_det_mod_pk_4block(Amod,pp,scratch) )
     {
      j=0;
      break;
     }
   }
  nmod_mat_clear(Amod);
  flint_free(scratch);
  if(j)
   j=_call_solve_Dixon(A);
  return j;
 }

#undef NDEBUG
#endif
