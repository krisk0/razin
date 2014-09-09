// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "nmod_mat_.h"

int 
nmod_mat_is_one(const nmod_mat_t m)
// m square. Return 1 iff m is identity
 {
  slong d=m->r;
  mp_limb_t mIJ;
  slong i,j;
  for(i=0;i<d;i++)
   {
    mp_limb_t* q=m->rows[i];
    for(j=0;j<d;j++)
     {
      mIJ = q[j];
      if( (i!=j) && mIJ )
       return 0;
      if( (i==j) && (mIJ!=1) )
       return 0;
     }
   }
  return 1;
 }
