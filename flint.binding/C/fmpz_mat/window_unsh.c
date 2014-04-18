// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// This file contains definition of fmpz_mat_window_unsh_struct and some methods
//  accessing it such as recursive triangular inverse computation

// fmpz_mat_window_unsh_struct is like tmod_mat_window_unsh_struct, but it is
//  for fmpz_mat rather than tmod_mat

typedef struct
 {
  slong r;
  slong c;
  fmpz** rows;
  slong delta;
 }
fmpz_mat_window_unsh_struct;
typedef fmpz_mat_window_unsh_struct fmpz_mat_window_unsh_t[1];
#define fmpz_mat_window_entry(mat,i,j) ((mat)->rows[(i)][(j)])

#include "flint/fmpz_mat.h"
mp_limb_t fmpz_to_t(const fmpz_t f);

static __inline__ void
fmpz_vec_init_set(fmpz* vec1, const fmpz* vec2, slong len2)
// ignore contents of vec1, perform deep copy
 {
  for ( ; len2--; )
   fmpz_init_set(vec1 + len2, vec2 + len2);
 }

static __inline__
void fmpz_clear_only(fmpz_t f)
 {
  if (COEFF_IS_MPZ(*f))
   _fmpz_clear_mpz(*f);
 }

static __inline__ void
fmpz_mat_to_window_unsh( fmpz_mat_window_unsh_t W, const fmpz_mat_t A )
// deep copy-constructor
 {
  slong i, cc=A->c, rc=A->r;
  W->delta = W->c = cc;
  W->r = rc;
  W->rows = flint_malloc( rc * sizeof(fmpz *) );
  fmpz* e = flint_malloc( rc * cc * sizeof(fmpz) );
  for( i=0; i<rc; i++, e += cc )
   {
    W->rows[i] = e;
    fmpz_vec_init_set(e, A->rows[i], cc);
   }
 }

static __inline__ void
fmpz_vec_kill_set(fmpz* tgt, fmpz* sou, slong len)
/*
destroy entries of tgt
shallow-copy entries of sou
*/
 {
  for( ; len--; )
   {
    fmpz_clear_only(tgt);
    *tgt++ = *sou++;
   }
 }

static __inline__ void
fmpz_mat_from_window_unsh( fmpz_mat_t R, fmpz_mat_window_unsh_t W )
/*
R dimensions should exactly match W dimensions
makes R virgin, copies all data from W to R, destroys W
do not call fmpz_mat_window_unsh_*clear() on W, after calling this procedure
*/
 {
  slong i, cc=W->c, rc=W->r, wd=W->delta;
  fmpz* T = R->entries;
  fmpz* S = W->rows[0];
  for (i = 0; i<rc; i++, S += wd)
   {
    R->rows[i]=T;
    fmpz_vec_kill_set( T, S, cc );
    T += cc;
   }
  flint_free( W->rows[0] );
  flint_free( W->rows );
 }

static __inline__ void
fmpz_mat_window_unsh_deep_clear( fmpz_mat_window_unsh_t W )
 {
  slong i, siz=W->r * W->c;
  fmpz* e = W->rows[0];
  for (i = 0; i < siz; i++, e++)
   fmpz_clear( e );
  flint_free( W->rows[0] );
  flint_free( W->rows );
 }

static __inline__ void
fmpz_mat_window_unsh_clear( fmpz_mat_window_unsh_t W )
 {
  flint_free( W->rows );
 }

static void
fmpz_mat_window_unsh_init(fmpz_mat_window_unsh_t W,
  const fmpz_mat_window_unsh_t S, slong r1, slong c1, slong r2, slong c2)
 {
  W->delta=S->delta;
  slong i, rD=r2-r1;
  W->rows = flint_malloc( rD * sizeof(fmpz*) );
  for (i = 0; i < rD; i++)
   W->rows[i] = S->rows[r1 + i] + c1;
  W->r = rD;
  W->c = c2 - c1;
 }

static __inline__ void
fmpz_clear_and_zero(fmpz_t f)
// works like fmpz_zero(), but does not assign zero twice --- one time
//  is quite enough
 {
  if (COEFF_IS_MPZ(*f))
   _fmpz_clear_mpz(*f);
  (*f) = WORD(0);
 }

static __inline__ void
fmpz_mat_window_unsh_zero( fmpz_mat_window_unsh_t W )
 {
  slong i, rc=W->r, cc=W->c, de=W->delta, j;
  fmpz* e = W->rows[0];
  for( i=0; i<rc; i++,e += de)
   for( j=0; j<cc; j++ )
    fmpz_clear_and_zero( e+j );
 }

#define MATR fmpz_mat_window_unsh_t
#define MUL_TYPE_ARG MATR R, const MATR A, const MATR B

static void fmpz_mat_window_unsh_print(const char* m,const MATR A)
 {
  #if 0
   printf("\n%s\n",m);
   slong i,j;
   for(i=0;i<A->r;i++)
    {
     for(j=0;j<A->c;j++)
      {
       //printf("%ld ",(slong)(fmpz_mat_window_entry(A,i,j)));
       fmpz_print( A->rows[i]+j ); printf(" ");
      }
     printf("\n");
    }
   printf("\n");
  #endif
 }

void
fmpz_neg_1arg( fmpz_t v )
// inspired by FLINT fmpz_neg
 {
  if (!COEFF_IS_MPZ(*v))
   *v = -*v;
  else
   {
    __mpz_struct* mpz_ptr = _fmpz_promote(v);
    // TODO: write effective subroutine to invert GMP mpz in-place
    mpz_neg(mpz_ptr, mpz_ptr);
   }
 }

mp_limb_t
fmpz_triU_inverse_dim2( MATR W )
 {
  fmpz* e0 = W->rows[0];
  fmpz* e1 = W->rows[1];
  mp_limb_t d = fmpz_to_t( e0 ) * fmpz_to_t( e1+1 );
  e1[0]=e1[1]; e1[1]=e0[0]; e0[0]=e1[0];
  e1[0]=WORD(0);
  fmpz_neg_1arg( e0+1 );
  return d;
 }

mp_limb_t
fmpz_triU_inverse_dim3( MATR W )
 {
  fmpz_mat_window_unsh_print("dim3 starts",W);
  fmpz* e0 = W->rows[0];
  fmpz* e1 = W->rows[1];
  fmpz* e2 = W->rows[2];
  mp_limb_t a0 = fmpz_to_t( e0 );
  mp_limb_t b0_b2 = fmpz_to_t( e1+1 ) * fmpz_to_t( e2+2 );
  fmpz_set_ui( e0, b0_b2 );
  fmpz_mat_window_unsh_print("UL fixed",W);
  // swap b2 <-> b0, using e2[0] as temp, don't zero e2[0]
  e2[0]=e2[2]; e2[2]=e1[1]; e1[1]=e2[0];
  fmpz* p=e1+2;
  fmpz_neg_1arg( p );
  fmpz_mat_window_unsh_print("B' counted",W);
  // e2[0] := a1
  e2[0]=e0[1];
  // e0[1] := vector (a1,a2) by 0th col of matrice B'
  fmpz* q=e0+1;
  fmpz_mul( q, e2, e1+1 );
  fmpz_neg_1arg( q );
  fmpz_mat_window_unsh_print("0th col of B' multiplied",W);
  // a2 := vector (a1,a2) by 1st col of matrice B'
  q++;p=e2+2;
  fmpz_mul   ( q,  q,  p );
  fmpz_addmul( q, e2, e1+2 );
  fmpz_neg_1arg( q );
  fmpz_mat_window_unsh_print("1st col of B' multiplied",W);
  // e2[0], thank you for being scratch
  e2[0]=WORD(0);
  // multiply B' by a0
  fmpz_mul_ui( p, p, a0 );
  p=e1+2;
  fmpz_mul_ui( p, p, a0 );
  p--;
  fmpz_mul_ui( p, p, a0 );
  fmpz_mat_window_unsh_print("B' multiplied by a0",W);
  return a0*b0_b2;
 }

static __inline__ void
fmpz_mat_window_unsh_mul_general_triu( MUL_TYPE_ARG )
// R := A*B where B is upper-triangular square
 {
  slong d0=R->r, d1=B->c, i,j,k;
  fmpz* p,* q,* rho;
  fmpz** r=R->rows;
  fmpz** b=B->rows;
  fmpz** a=A->rows;
  slong B_delta=B->delta;
  fmpz_t B_00; fmpz_init_set( B_00, b[0] );
  fmpz_t scr; fmpz_init( scr );
  for( i=0; i<d0; i++ )
   {
    rho=r[i];
    p=a[i];
    fmpz_mul( rho, p, B_00 );
    rho++;
    for( j=1; j<d1; j++ )
     {
      fmpz_clear_and_zero( rho );
      q=b[0]+j;
      // multiply j entries
      for( k=0; k<j; k++,q += B_delta  )
       fmpz_addmul( rho, p+k, q );
      // multiply j'th entry
      fmpz_addmul( rho, p+j, q );
      rho++;
     }
   }
  fmpz_clear( scr );
  fmpz_clear( B_00 );
 }

static __inline__ void
fmpz_mat_window_unsh_mul_negate_triu_general( MUL_TYPE_ARG )
// R := -A*B where A is upper-triangular square
 {
  slong d0=R->r, d1=B->c, i,j,k;
  slong d0_minus=d0-1;
  fmpz* p,* q,* rho;
  fmpz** r=R->rows;
  fmpz** b=B->rows;
  fmpz** a=A->rows;
  slong B_delta=B->delta;
  fmpz_t scr; fmpz_init( scr );
  // All rows but last
  for (i=0; i<d0_minus; i++)
   {
    rho=r[i];
    p=a[i];     // to skip i zeroes at start of a, add i
    for ( j=0; j<d1; j++ )
     {
      q=b[i]+j; // skip i rows of b, then go to j'th column
      fmpz_clear_and_zero( scr );
      for( k=i; k<d0; k++, q += B_delta )
       fmpz_addmul( scr, p+k, q );
      fmpz_neg_1arg( scr );
      fmpz_swap( rho, scr ); // if R entry was small, scr becomes small
      rho++;
     }
   }
  // last row: scalar multiply
  p=r[d0_minus]; q=b[d0_minus];
  fmpz_set( scr, a[d0_minus]+d0_minus ); fmpz_neg_1arg( scr );
  for( i=d1; i--; )
   {
    fmpz_mul( p, q, scr );
    p++;
    q++;
   }
  fmpz_clear( scr );
 }

static __inline__ void
fmpz_mat_triU_C_transform( MATR R, slong A_dim, slong D_dim,
  const MATR A, MATR Cscr, const MATR D)
// TODO: Use asymptotically fast multiplication  
 {
  MATR C;
  fmpz_mat_window_unsh_init( C, R, 0, A_dim, A_dim, A_dim+D_dim );
  fmpz_mat_window_unsh_print("A=",A);
  fmpz_mat_window_unsh_print("C=",C);
  fmpz_mat_window_unsh_mul_negate_triu_general( Cscr, A, C ); // Cscr = -A*C
  fmpz_mat_window_unsh_print("-A*C=",Cscr);
  fmpz_mat_window_unsh_mul_general_triu       ( C, Cscr, D ); // C = Cscr * D
  fmpz_mat_window_unsh_print("D=",D);
  fmpz_mat_window_unsh_print("*D=",C);
  fmpz_mat_window_unsh_clear(C);
 }

static __inline__ void
fmpz_mat_window_unsh_mult_ui( MATR R, mp_limb_t m )
// R := R*m
 {
  slong rc=R->r, cc=R->c, i,j;
  fmpz* p;
  for( i=rc; i--; )
   {
    p=R->rows[i];
    for( j=cc; j--; p++ )
     fmpz_mul_ui( p, p, m );
   }
 }

mp_limb_t
fmpz_triU_inverse_rec_smallDet(MATR alpha)
 {
  slong m=alpha->r;
  assert(m>1);
  if(m==2)
   return fmpz_triU_inverse_dim2( alpha );
  if(m==3)
   return fmpz_triU_inverse_dim3( alpha );
  /*
  Use recursive formula
 
  (A C)'    =   (  A'  -A' C D' )
  (  D)         (          D'   )
  */
  slong n0=m>>1;
  slong n1=m-n0;
  MATR A, S, D;
  fmpz_mat_window_unsh_init( A, alpha, 0,  0,  n0, n0 );
  fmpz_mat_window_unsh_init( D, alpha, n0, n0,  m,  m );
  mp_limb_t d0=fmpz_triU_inverse_rec_smallDet( A );
  mp_limb_t d1=fmpz_triU_inverse_rec_smallDet( D );
  // need temp storage to replace C with product A'*C*D', using zeroes at
  //  lower-left corner
  fmpz_mat_window_unsh_init( S, alpha, n1, 0,   m, n1 );
  fmpz_mat_triU_C_transform( alpha, n0, n1, A, S, D );
  // clean rubber in lower-left corner, forget it
  fmpz_mat_window_unsh_zero(S);         fmpz_mat_window_unsh_clear(S);
  fmpz_mat_window_unsh_mult_ui( A, d1 );
  fmpz_mat_window_unsh_mult_ui( D, d0 );
  fmpz_mat_window_unsh_clear(D);
  fmpz_mat_window_unsh_clear(A);
  return d0*d1;
 }

#undef MUL_TYPE_ARG
#undef MATR

void
fmpz_triU_inverse_smallDet(fmpz_mat_t T, fmpz_t d, const fmpz_mat_t S)
/*
S: square upper-triangular, count of rows >= 4, diagonal positive
det S is in range 2 .. 2**64-1
upon exit T,d are such that T/d*source S=identity, d divides det S
aliasing allowed
*/
 {
  mp_limb_t de;
  fmpz_mat_window_unsh_t W;
  fmpz_mat_to_window_unsh( W, S );
  de=fmpz_triU_inverse_rec_smallDet( W );
  fmpz_mat_from_window_unsh( T, W );
  // W has been eaten, no need to un-allocate
  fmpz_set_ui( d, de );
 }
