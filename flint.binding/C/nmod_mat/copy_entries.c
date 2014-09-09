// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "nmod_mat_.h"

void 
nmod_mat_copy_entries(nmod_mat_t r,const nmod_mat_t s)
 {
  slong max_i=s->r,i;
  size_t size=s->c*sizeof(mp_limb_t);
  for(i=0;i<max_i;i++)
   memcpy( r->rows[i], s->rows[i], size );
 }
