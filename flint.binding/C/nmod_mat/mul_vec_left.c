// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include "../ulong_extras/ulong_extras_.h"
#include "../ulong_extras/longlong_.h"
#include "nmod_mat_.h"

void 
nmod_mat_mul_vec_left(mp_limb_t* r, const mp_limb_t* s, const nmod_mat_t m)
// r = s*A, not optimized
 {
  slong i_max=m->c,j_max=m->r,i,j;
  mp_limb_t* rP;
  for(i=0,rP=r;i<i_max;i++,rP++)
   {
    // r[i] = s * m column no. i
    VECTOR_DOT_HEAD_greedy( s[0], m->rows[0][i] );
    for(j=1;j<j_max;j++)
     VECTOR_DOT_BODY_greedy( s[j], m->rows[j][i] );
    #if SPEEDUP_NMOD_RED3
     VECTOR_DOT_TAIL( rP[0], m->mod.n,m->mod.ninv,m->mod.norm );
    #else
     VECTOR_DOT_TAIL_3arg( rP[0], m->mod.n,m->mod.ninv );
    #endif
   }
 }
