// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>

mp_limb_t nmod_mat_det_mod_pk_4block(nmod_mat_t M,const p_k_pk_t pp,mp_limb_t* scrtch);
void nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim);

mp_limb_t
fmpz_mat_det_mod_pk_4block(fmpz_mat_t M,nmod_mat_t Mmod,const p_k_pk_t pp,
  mp_limb_t* scrtch)
/*
M: square of dimension >= 1
return determinant of M modulo p**k
modify Mmod and scrtch
*/
 {
  _nmod_mat_set_mod(Mmod, pp.p_deg_k);
  Mmod->mod.n <<= Mmod->mod.norm;  // make sure mod.n has higher bit set
  #if SPEEDUP_NMOD_RED3
   mp_limb_t t= - Mmod->mod.n;
   Mmod->mod.norm = n_mulmod_preinv_4arg( t,t, Mmod->mod.n,Mmod->mod.ninv );
  #else
   Mmod->mod.norm=0;
  #endif
  fmpz_mat_get_nmod_mat(Mmod, M);
  return nmod_mat_det_mod_pk_4block(Mmod,pp,scrtch);
 }

mp_limb_t
fmpz_mat_det_mod_pk_3arg(fmpz_mat_t M,mp_limb_t p,ulong k)
 {
  slong dim=M->c;
  p_k_pk_t pp; pp.p=p; pp.k=k; pp.p_deg_k=n_pow(p,k);
  mp_limb_t* scratch;
  if(dim>4)
   scratch=flint_malloc( 4*(dim-4)*sizeof(mp_limb_t) );
  nmod_mat_t Mmod;
  nmod_mat_init_square_2arg(Mmod, dim);
  mp_limb_t r=fmpz_mat_det_mod_pk_4block(M,Mmod,pp,scratch);
  nmod_mat_clear(Mmod);
  if(dim>4)
   flint_free(scratch);
  return r;
 }
