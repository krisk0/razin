// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// No, code below is not obfuscated.
// When I write a program that needs to run fast, it comes this way

// This file is auto-generated from C/nmod_mat/nmod_mat_HNF-debug.c

#define NDEBUG 1

#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

static __inline__ slong 
DKryskov_nmod_find_nonzero(nmod_mat_t A,slong col,mp_limb_t det_tgt)
/*
 -3 if diagonal element is 1
 -2 if all is zero
 -1 if single non-zero
 line>0 if A[line,i] is non-zero
TODO: multiply by -1 to get smaller element
*/
 {
  slong i0=-1,i1=-1,k=col,max_k=A->r;
  mp_limb_t m0=det_tgt;
  mp_limb_t m1=m0;
  mp_limb_t* kP;
  mp_limb_t Ak;
  while(k<max_k)
   {
    kP = &nmod_mat_entry(A, k, col);
    Ak = *kP;
    if(Ak == 0)
     {
      k += 1;
      continue;
     }
    if(Ak > det_tgt) 
     {
      Ak %= det_tgt;
      *kP=Ak;
      if(Ak == 0)
       {
        k += 1;
        continue;
       }
     }
    // new champion?
    if(Ak < m0)
     {
      m1=m0; i1=i0;
      m0=Ak;  i0=k;
      if( (m0 == 1) && (i1 != -1) )
       {
        if( k > col )
         MP_PTR_SWAP( A->rows[col], A->rows[i0] );
        return -3;
       }
      k += 1;
      continue;
     }
    // new silver pretender?
    if(Ak < m1)
     {
      m1=Ak; i1=k;
      if( m0 == 1 )
       {
        if( i0 > col )
         MP_PTR_SWAP( A->rows[col], A->rows[i0] );
        return -3;
       }
     }
    k += 1;
   }
  if(i1 > -1)
   {
    if( (i0 != col) && (i1 != col) )
     {
      // None of the winners at the pedestal, fix it
      MP_PTR_SWAP( A->rows[col], A->rows[i0] );
      return i1;
     }
    // One of the winners already here, refer to another
    if(i0 == col)
     return i1;
    return i0;
   }
  // Not more than one non-zero
  if(i0 == -1)
   return -2;
  if(i0 != col)
   MP_PTR_SWAP( A->rows[col], A->rows[i0] );
  return -1;
 }

#include <assert.h>
static __inline__ mp_limb_t
DKryskov_gcd_ui(mp_limb_t* u,mp_limb_t* v,mp_limb_t x,mp_limb_t y,mp_limb_t n)
//return g=gcd(x,y) and set numbers u,v such that u*x+v*y=g modulo n
 {
  mp_limb_t g;
  if(x<y)
   {
    // TODO: n_xgcd is in C and GMP mpn_gcd_1 is in ASM.
    // Should I use a GMP subroutine instead of n_xgcd?
    g=n_xgcd(v,u,y,x);
    assert( g == (*v)*y - (*u)*x );
    *u = n-( (*u) % n );
    *v %= n;
    return g;
   }
  g=n_xgcd(u,v,x,y);
  assert( g == (*u)*x - (*v)*y );
  *v = n-( (*v) % n );
  *u %= n;
  return g;
 }

static __inline__ slong
DKryskov_nmod_zero_line(nmod_mat_t A,slong i,slong j,mp_limb_t n,
 mp_limb_t* scrth)
//returns 1 iff new A[i,i] becomes 1
 {
  assert(i != j);
  slong iPLUS=i+1;
  mp_limb_t* alpha=A->rows[i]+iPLUS;
  mp_limb_t* betta=A->rows[j]+iPLUS;
  mp_limb_t x=alpha[-1];
  mp_limb_t y=betta[-1];
  assert(x>1);
  assert(y>1);
  // TODO: multiply 2nd line by -1 so gcd runs faster
  mp_limb_t u,v,g;
  g=DKryskov_gcd_ui(&u,&v,x,y,n);
  assert(g);
  mp_limb_t vec_len=A->c-iPLUS;
  _nmod_vec_scalar_mul_nmod(    scrth, alpha, vec_len, u  , A->mod );
  _nmod_vec_scalar_addmul_nmod( scrth, betta, vec_len, v   , A->mod );
  _nmod_vec_scalar_mul_nmod(    betta, betta, vec_len, (x/g) , A->mod );
  _nmod_vec_scalar_addmul_nmod( betta, alpha, vec_len, (n-y/g) , A->mod );
  memcpy(alpha,scrth,vec_len*sizeof(mp_limb_t));
  alpha[-1]=g;
  betta[-1]=0;
  return g==1;
 }

static __inline__ void 
DKryskov_nmod_Gauss_upper_last_col(nmod_mat_t A,slong last_col)
 {
  mp_limb_t s=nmod_mat_entry( A, last_col, last_col );
  assert(s);
  slong j;
  for(j=last_col;j--;)
   nmod_mat_entry( A, j, last_col ) %= s;
 }

static __inline__ void 
DKryskov_nmod_Gauss_upper(nmod_mat_t A)
/*
zap A over diagonal
attempt to avoid row operations if possible
*/
 {
  slong last_i=A->c-1;
  slong j,i;
  mp_limb_t n=nmod_mat_entry(A,1,1);
  for(i=2;i<=last_i;i++)
   n *= nmod_mat_entry(A,i,i);
  mp_limb_t* sP,*tP;
  mp_limb_t s,t_ori,t_upd;
  for(i=1;i<last_i;i++)
   {
    sP=&nmod_mat_entry( A, i,i );
    mp_limb_t s=*sP;
    assert(s);
    slong v_len=A->c-i;
    for(j=i;j--;)
     {
      tP=&nmod_mat_entry( A, j,i );
      t_ori=*tP;
      if(t_ori >= s)
       {
        t_upd=t_ori % n;
        if(t_upd < s)
         // don't need vector operation
         *tP = t_upd;
        else
         {
          mp_limb_t q=t_upd/s;
          _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-q, A->mod );
          *tP %= n;
         }
       }
      assert( *tP < s );
     }
    n /= s; 
   }
  DKryskov_nmod_Gauss_upper_last_col( A, last_i );
 }

static __inline__ void 
DKryskov_nmod_early_abort(nmod_mat_t A,slong e)
//all columns after e have 1 on diagonal, so save some effort
// TODO: apply decreasing modulo to avoid row operations
 {
  slong m=A->c;
  slong i,j;
  mp_limb_t* tP,* sP;
  slong ePLUS=e+1;
  j=sizeof(mp_limb_t)*(m-ePLUS);
  for(i=m;i--;)
   {
    tP=A->rows[i]+ePLUS;
    memset(tP,0,j);
    if(i>=ePLUS)
     tP[i-ePLUS]=1;
   }
  // now like Gauss_upper, but don't change anything after column e
  mp_limb_t n=A->mod.n;
  for(i=1;i<e;i++)
   { 
    sP=&nmod_mat_entry( A, i, i );
    mp_limb_t s=*sP;
    assert(s);
    slong v_len=ePLUS-i;
    for(j=i;j--;)
     {
      tP=&nmod_mat_entry( A, j, i );
      mp_limb_t  t=*tP;
      if(t)
       {
        mp_limb_t q=t/s;
        if(q)
         _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-q, A->mod );
        assert( *tP < s);
       }
     }
   }
  DKryskov_nmod_Gauss_upper_last_col( A, e );
 }

static __inline__ void
DKryskov_nmod_reduce_diag(nmod_mat_t A,slong i,mp_limb_t det_tgt,mp_limb_t* scratch)
 {
  assert(i<A->c-1);
# if 1
   if( det_tgt % nmod_mat_entry(A,i,i) )
    {
     assert( 0 == nmod_mat_entry(A,i+1,i) % det_tgt );
     nmod_mat_entry(A,i+1,i)=det_tgt;
     (void)DKryskov_nmod_zero_line(A,i,i+1,det_tgt,scratch);
    }
# else
   // code below is incorrect: bad algorithm
   mp_limb_t aII=nmod_mat_entry(A,i,i);
   if( det_tgt % aII )
    {
     nmod_mat_entry(A,i,i)=n_gcd_full(det_tgt,aII);
    }
# endif
 }

static __inline__ void
DKryskov_nmod_reduce_last( mp_limb_t* se_corner, mp_limb_t det_tgt )
 {
  *se_corner = n_gcd_full( *se_corner, det_tgt );
 }

static __inline__ void
DKryskov_nmod_1_lower(nmod_mat_t A,slong col,slong j,mp_limb_t n)
/*
A[col,col] is known to be 1
zap column col starting from row j
operate modulo n
*/
 {
  slong m=A->c;
  mp_limb_t* sP=A->rows[col]+col;
  assert( 1 == *sP );
  //TODO: skip counting elements of column col, to save time
  slong v_len = m-col;
  while(j < m)
   {
    mp_limb_t* tP=A->rows[j]+col;
    mp_limb_t t_ori=*tP;
    if(t_ori)
     {
      mp_limb_t t_upd=t_ori % n;
      if(t_upd)
       {
        _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-t_upd, A->mod );
       }
      else
       *tP = 0; // TODO: this line can be deleted?
     }
    ++j;
   }
 }

slong 
nmod_mat_HNF(nmod_mat_t A)
/*
HNF: Hermite Normal Form as defined by W.Stein/C.Pernet

Input data:
~~~~~~~~~~
A: square, non-singular over integers, d=det(A) counted over Z,
 n=A->mod.n is a multiple of abs(d) (n >= abs(d), abs(d) divides n)

Output data:
~~~~~~~~~~~
Define H := HNF( input A )

Output A contains matrice H: its diagonal and above entries match those of H:

for i<=j A[i,j]=H[i,j]
 
Return code: 0
 
Literature:
~~~~~~~~~~
"Hermite Normal Form Computation Using Modulo Determinant Arithmetic"
aka DomichKannanTrotter87.pdf 
*/
 {
  slong m=A->c;
  slong i;
  slong bad_n=0;
  mp_limb_t* scratch=(mp_limb_t*)malloc( sizeof(mp_limb_t) * (m-1) );
  mp_limb_t det_tgt=A->mod.n;
  for(i=0;;i++)
   {// main loop: zap lower, fix diagonal, maintain decreasing modulo
    assert( (i>=0) && (i<m) );
    slong j=DKryskov_nmod_find_nonzero(A,i,det_tgt);
    if(j==-2)
     {
      nmod_mat_entry(A,i,i)=det_tgt;
      if(i<m-1)
       {
        free(scratch);
        DKryskov_nmod_early_abort(A,i);
        return 0;
       }
      // success, diagonal is correct
      break;
     }
    if(j==-1)
     {
      if( det_tgt % nmod_mat_entry(A,i,i) )
       {
        if(i==m-1)
         {
          // bad input: n on input not a multiple of matrice determinant
          DKryskov_nmod_reduce_last( &nmod_mat_entry(A,i,i), det_tgt );
          bad_n=1;
          break;
         }
       }
      if( i==m-1 )
       //skip to Gauss_upper
       break; 
     }
    if(j==-3)
     { 
      // some rows operations required, diagonal element = 1
      DKryskov_nmod_1_lower(A,i,i+1,det_tgt);
      // no need to fix diagonal or change det_tgt
      continue;
     }
    if(j>=0)
     {
      // some rows operations required, general case
      if(DKryskov_nmod_zero_line(A,i,j,det_tgt,scratch))
       {
        j=i+1;
        if(j<m)
         DKryskov_nmod_1_lower(A,i,j,det_tgt);
        continue;
       }
      j=i;
      while( ++j < m)
       {
        if( nmod_mat_entry(A,j,i) % det_tgt )
         {
          if(DKryskov_nmod_zero_line(A,i,j,det_tgt,scratch))
           {
            if( ++j < m )
             DKryskov_nmod_1_lower(A,i,j,det_tgt);
            break;
           }
         }
       }
     }
    DKryskov_nmod_reduce_diag(A,i,det_tgt,scratch);
    assert( 0 == det_tgt % nmod_mat_entry(A,i,i) );
    det_tgt /= nmod_mat_entry(A,i,i);
    if( (det_tgt==1) && (i<m-1) )
     {
      free(scratch);
      DKryskov_nmod_early_abort(A,i);
      return 0;
     }
   }// main loop end
  free(scratch);
  DKryskov_nmod_Gauss_upper(A);
  return bad_n;
 }

#include "C/nmod_mat/nmod_mat_HNF_nonsq.c"

slong
nmod_mat_HNF_nonsquare(nmod_mat_t A)
/*
Suppose matrice B exists such that n=A=>mod.n>1, B mod n=A, c=A->c <= A->r, 
 rank of B over Z equals c, H=upper c rows of HNF(B), n divides det(H)>0
 
compute H'=upper c rows of HNF of matrice B modulo a lattice with determinant 
 n, chosen in such a way that det H'=n
        
Result H' returned in diagonal and over-diagonal entries of A
 
Return code:0
*/
 {
  slong c=A->c;
  slong i,j;
  mp_limb_t* scratch=(mp_limb_t*)malloc( sizeof(mp_limb_t) * (c-1) );
  mp_limb_t det_tgt=A->mod.n;
  for(i=0;;i++)
   {// main loop: zap lower, fix diagonal, maintain decreasing modulo
    assert( (i>=0) && (i<c) );
    j=DKryskov_nmod_find_nonzero(A,i,det_tgt);
    if(i == c-1)
     {
      free(scratch);
      DKryskov_HNF_zap_below_last_col(A,i,j,det_tgt);
      j=DKryskov_fix_diagonal_tail(A,i,det_tgt);
      DKryskov_nmod_Gauss_upper(A);
      return j;
     }
    DKryskov_HNF_zap_below(A,i,j,det_tgt,scratch);
    DKryskov_nmod_reduce_diag(A,i,det_tgt,scratch);
    det_tgt /= nmod_mat_entry(A,i,i);
   }
 }


#undef NDEBUG

/*
I would be greatly irritated if you tell me that the algorithms implemented in
 this file are wrong WITHOUT GIVING EXAMPLE OF INPUT DATA that make them fail

Report bugs via Github mechanism or e-mail

My e-mail is in my blog, detailed information on how to get it is close to tail
 of setup.py
*/