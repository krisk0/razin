// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include "../tmod_mat/tmod_mat_.h"

mp_limb_t fmpz_to_t(const fmpz_t);

void nmod_mat_mod_t(nmod_mat_t tgt, const fmpz_mat_t sou)
/*
 initialize tgt like nmod_mat_init(), fill mod with zeroes, 
  compute tgt := sou modulo 2**64
*/
 {
  slong rc=sou->r, cc=sou->c, i, j;
  tgt->r=rc; tgt->c=cc;
  nmod_t* mod=&(tgt->mod);
  mod->n=mod->norm=mod->ninv=0;
  mp_limb_t* e = (mp_limb_t*)flint_malloc( rc * cc * sizeof(mp_limb_t) );
  tgt->entries=e;
  tgt->rows = (mp_limb_t**)flint_malloc( rc * sizeof(mp_limb_t*) );
  fmpz* s;
  for(i=0;i<rc;i++,e += cc)
   {
    tgt->rows[i]=e; 
    s=sou->rows[i];
    for(j=0;j<cc;j++)
     e[j] = fmpz_to_t( s+j );
   }
 }
