// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

void invert_permutation(mp_limb_t* r,mp_limb_t* s,slong n)
 {
  slong i;
  for(i=0;i<n;i++)
   r[s[i]] = i;
 }
