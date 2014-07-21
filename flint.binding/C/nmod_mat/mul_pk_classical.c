// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include "../ulong_extras/ulong_extras_.h"
#include "../ulong_extras/longlong_.h"
#include "nmod_mat_.h"

void
nmod_mat_mul_pk_classical(nmod_mat_t R,nmod_mat_t A,nmod_mat_t B)
 {
  slong i_max=R->r,j_max=R->c,k_max=B->r,i,j,k;
  mp_limb_t** const Rrows=R->rows;
  mp_limb_t** const Arows=A->rows;
  mp_limb_t** const Brows=B->rows;
  mp_limb_t* rho;
  mp_limb_t* alpha;
  mp_limb_t* betta;
  const mp_limb_t n=R->mod.n;
  const mp_limb_t ninv=R->mod.ninv;
  #if SPEEDUP_NMOD_RED3
   const mp_limb_t two_128_mod_n=R->mod.norm;
  #endif
  for(i=0;i<i_max;i++)
   {
    rho=Rrows[i];
    alpha=Arows[i];
    for(j=0;j<j_max;j++)
     {
      VECTOR_DOT_HEAD( alpha[0], Brows[0][j] );
      for(k=1;k<k_max;k++)
       VECTOR_DOT_BODY( alpha[k], Brows[k][j] );
      #if SPEEDUP_NMOD_RED3
       VECTOR_DOT_TAIL( rho[j], n,ninv,two_128_mod_n );
      #else
                                               // for 4x4x100 multiplication,
       VECTOR_DOT_TAIL_3arg( rho[j], n,ninv ); //  17% slower
      #endif
     }
   }
 }
