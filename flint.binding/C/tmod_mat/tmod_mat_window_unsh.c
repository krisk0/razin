// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// some code borrowed from FLINT nmod_mat/*.c, however data layout is different.
//  Notice delta field

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

static void
tmod_mat_to_window_unsh(tmod_mat_window_unsh_t W, const tmod_mat_t A,
 slong r1, slong c1, slong r2, slong c2)
 {
  W->delta=A->c;
  slong i, rD=r2-r1;
  W->rows = flint_malloc( rD * sizeof(mp_limb_t *) );
  for (i = 0; i < rD; i++)
   W->rows[i] = A->rows[r1 + i] + c1;
  W->r = rD;
  W->c = c2 - c1;
 }

static void 
tmod_mat_window_unsh_clear( tmod_mat_window_unsh_t W )
 {
  flint_free( W->rows );
 }

static __inline__ void 
tmod_mat_invert_tril_dim3( tmod_mat_window_unsh_t W )
/*
            1                                      1
inverse of w0  1     according to sympy equals   -w0       1
           w1 w2  1                             w0*w2-w1  -w2    1
*/
 {
  // input data
  mp_limb_t* e1=W->rows[1],* e2=W->rows[2];
  mp_limb_t w0=e1[0], w1=e2[0], w2=e2[1];
  // calculate and output
  e1[0] = -w0;
  e2[0] = w0*w2-w1;  e2[1] = -w2;
 }

static void
tmod_mat_invert_tril( tmod_mat_window_unsh_t alpha )
/*
on entry alpha: square, lower-triangular, diagonal 1
on exit  alpha: old alpha inverted
*/
 {
  long m=alpha->r;
  if(m==1)
   return;
  if(m==2)
   {
    mp_limb_t* e=alpha->rows[1];
    *e = -*e;
    return;
   }
  if(m==3)
   {
    tmod_mat_invert_tril_dim3( alpha );
    return;
   }
  exit(1);
 }
