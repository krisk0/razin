// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/nmod_mat.h>

void
nmod_mat_transpose_square_tgt_virgin(nmod_mat_t tgt,const nmod_mat_t sou)
 {
  slong n=tgt->r,i,j;
  const mp_limb_t** const sou_rows=sou->rows;
  mp_limb_t* p;
  mp_limb_t* q;
  for(i=0;i<n;i++)
   {
    p=sou_rows[i];
    q=tgt->entries+i;
    for(j=0;j<n;j++,p++,q += n)
     *q = *p;
   }
 }
