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
 p += delta;
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
tmod_mat_window_unsh_init(tmod_mat_window_unsh_t W,
  const tmod_mat_window_unsh_t S, slong r1, slong c1, slong r2, slong c2)
 {
  W->delta=S->delta;
  slong i, rD=r2-r1;
  W->rows = flint_malloc( rD * sizeof(mp_limb_t *) );
  for (i = 0; i < rD; i++)
   W->rows[i] = S->rows[r1 + i] + c1;
  W->r = rD;
  W->c = c2 - c1;
 }

static void
tmod_mat_window_unsh_mirror( tmod_mat_window_unsh_t M, slong rc, slong cc )
// allocate matrice and pretend that it is tmod_mat_window_unsh_struct
 {
  slong i;
  M->rows = flint_malloc( rc * sizeof(mp_limb_t*) );
  mp_limb_t* e = M->rows[0] = flint_malloc(rc * cc * sizeof(mp_limb_t));
  M->r = rc;
  M->delta = M->c = cc;
  for ( i=0; i<rc; i++, e += cc )
    M->rows[i] = e;
 }

static void 
tmod_mat_window_unsh_clear( tmod_mat_window_unsh_t W )
 {
  flint_free( W->rows );
 }

static void 
tmod_mat_window_unsh_unmirror( tmod_mat_window_unsh_t M )
 {
  flint_free( M->rows[0] );
  flint_free( M->rows );
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

#define MATR tmod_mat_window_unsh_t
#define MUL_TYPE_ARG MATR R, const MATR A, const MATR B

static __inline__ void
tmod_mat_window_unsh_mul_negate_tril1_general( MUL_TYPE_ARG )
// Count R=-A*B where A is square lower-tria with 1 on diagonal
 {
  slong d0=R->r, d1=B->c, i,j,k;
  mp_limb_t* p,* q,* rho;
  mp_limb_t** r=R->rows;
  mp_limb_t** b=B->rows;
  mp_limb_t** a=A->rows;
  mp_limb_t sum;
  slong B_delta=B->delta;
  // 0th row copy-negate
  p=r[0]; q=b[0];
  for (i=d1; i--; )
   *p++ = -*q++;
  // other rows
  for (i=1; i<d0; i++)
   {
    rho=r[i];
    for ( j=0; j<d1; j++ )
     {
      p=a[i]; 
      q=b[0]+j; 
      sum=0;
      for ( k=0; k<i; k++, q += B_delta )
       sum += p[k] * q[0];
      sum += q[0];   // multiply by 1
      *rho++ = -sum; // negate, store, shift rho by 1
     }
   }
 }

static __inline__ void
tmod_mat_window_unsh_mul_general_tril1( MUL_TYPE_ARG )
// Count R=A*B where B is square lower-tria with 1 on diagonal
 {
  slong d0=A->r, d1=B->c, i,j,k,j_plus;
  mp_limb_t* p,* q,* rho;
  slong d1_minus=d1-1;
  mp_limb_t** r=R->rows;
  mp_limb_t** b=B->rows;
  mp_limb_t** a=A->rows;
  slong R_delta=R->delta;
  slong A_delta=A->delta;
  slong B_delta=B->delta;
  mp_limb_t sum;
  for (i=0; i<d0; i++)
   {
    rho=r[i];  // will form i-th row of result, all columns but last
    for (j=0; j<d1_minus; j++ )
     {
      j_plus=j+1;
      p = a[i]+j_plus;
      q = b[0]+j+j_plus*B_delta; // skip 1 or more upper entries of j-th column
      sum = p[-1];              // multiplication by 1
      for (k=d1-j_plus; k--; p++, q += B_delta)
       sum += p[0] * q[0];
      *rho++ = sum;
     }
   }
  // copy last column
  p = r[0]+d1_minus; q=a[0]+d1_minus;
  for( i=d0; i--; p += R_delta, q += A_delta )
   *p = *q;
 }

static __inline__ void
tmod_mat_invert_tril_C_transform( MATR R, slong A_dim, slong D_dim,
  const MATR A, MATR Cscr, const MATR D)
 {
  MATR C;
  tmod_mat_window_unsh_init( C, R, A_dim, 0, A_dim+D_dim, A_dim );
  tmod_mat_window_unsh_mul_negate_tril1_general( Cscr, D, C ); // Cscr = -D*C
  tmod_mat_window_unsh_mul_general_tril1       ( C, Cscr, A ); // C = Cscr * A
  tmod_mat_window_unsh_clear(C);
 }

static void
tmod_mat_invert_tril( MATR alpha )
/*
on entry alpha: square, lower-triangular, diagonal of 1
on exit  alpha: old alpha inverted
*/
 {
  slong m=alpha->r;
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
  /*
  Use recursive formula
       
  (A  )'    =   (     A'          )
  (C D)         ( -D' C A'     D' )
  */
  slong n0=m>>1;
  slong n1=m-n0;
  MATR A, C, D;
  tmod_mat_window_unsh_init( A, alpha,  0,  0, n0, n0 );
  tmod_mat_window_unsh_init( D, alpha, n0, n0,  m,  m );
  tmod_mat_invert_tril( A );
  tmod_mat_invert_tril( D );
  // need temp storage to replace C with product D'*C*A'
  tmod_mat_window_unsh_mirror( C, n1, n0 );
  tmod_mat_invert_tril_C_transform( alpha, n0, n1, A, C, D );
  tmod_mat_window_unsh_unmirror(C);
  tmod_mat_window_unsh_clear(D);
  tmod_mat_window_unsh_clear(A);
 }

#undef MUL_TYPE_ARG
#undef MATR
