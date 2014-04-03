// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

typedef struct
 {
  slong r;
  slong c;
  mp_limb_t** rows;
  slong delta;      
 }
tmod_mat_window_unsh_struct;
typedef tmod_mat_window_unsh_struct tmod_mat_window_unsh_t[1];
#define tmod_mat_window_entry(mat,i,j) ((mat)->rows[(i)][(j)])

#if 0
 // shift by delta to get next in column
 mp_limb_t* p=mat->rows[0]+col; 
 assert( p[0] == tmod_mat_window_entry(mat,0,col) );
 p[0] += delta;
 assert( p[0] == tmod_mat_window_entry(mat,1,col) );
#endif
