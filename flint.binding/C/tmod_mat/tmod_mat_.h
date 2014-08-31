// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef TMOD_MAT__H
#define TMOD_MAT__H

#include <flint/flint.h>

typedef struct
 {
  slong r;
  slong c;
  mp_limb_t**   rows;
  mp_limb_t*  entries;
 }
tmod_mat_struct;
typedef tmod_mat_struct tmod_mat_t[1];
#define tmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])

#endif
