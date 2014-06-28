// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
mp_limb_t nmod_mat_det_mod_pk(nmod_mat_t M,mp_limb_t p,ulong k,
 mp_limb_t p_deg_k,mp_limb_t* scrtch);
void nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim);

mp_limb_t
fmpz_mat_det_mod_pk(fmpz_mat_t M,nmod_mat_t Mmod,
           mp_limb_t p,mp_limb_t p_deg_k,ulong k,
           mp_limb_t* scrtch)
/*
M: square of dimension >= 1
return determinant of M modulo p**k
modify Mmod and scrtch
*/
 {
  _nmod_mat_set_mod(Mmod, p_deg_k);
  Mmod->mod.n <<= Mmod->mod.norm;  // make sure mod.n has higher bit set
  Mmod->mod.norm = 0;              // now normalizer is zero
  fmpz_mat_get_nmod_mat(Mmod, M);
  return nmod_mat_det_mod_pk(Mmod,p,k,p_deg_k,scrtch);
 }

mp_limb_t
fmpz_mat_det_mod_pk_3arg(fmpz_mat_t M,mp_limb_t p,ulong k)
 {
  slong dim=M->c;
  mp_limb_t pk=n_pow(p,k);
  mp_limb_t* scratch;
  if(dim>4)
   scratch=flint_malloc( 4*(dim-4)*sizeof(mp_limb_t) );
  nmod_mat_t Mmod;
  nmod_mat_init_square_2arg(Mmod, dim);
  mp_limb_t r=fmpz_mat_det_mod_pk(M,Mmod,p,pk,k,scratch);
  nmod_mat_clear(Mmod);
  if(dim>4)
   flint_free(scratch);
  return r;
 }
