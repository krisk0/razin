#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#define BUG0_nmod_mat_HNF 1

#if 0
 WARNING

  Two files are put online in same directory: 
   nmod_mat_HNF-debug.c (this file) and mod_mat_HNF.c

  However the latter is auto-generated out of the former with 
   ./delete_debug_C_code.sh
#endif

static __inline__ long 
DKryskov_nmod_find_nonzero(nmod_mat_t A,long col,mp_limb_t det_tgt)
 {
  long i0=-1,i1=-1,k=col;
  mp_limb_t m0=A->mod.n;
  mp_limb_t m1=m0;
  mp_limb_t* kP;
  mp_limb_t Ak;
  while(k<A->c)
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
      k += 1;
      continue;
     }
    // new silver pretender?
    if(Ak < m1)
     {
      m1=Ak; i1=k;
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
   {
    MP_PTR_SWAP( A->rows[col], A->rows[i0] );
   }
  return -1;
 }

static void nmod_mat_print(nmod_mat_t A)
// nmod_mat_print_pretty was printing garbage like [%%2lu %%2lu]
 {
  long i,j;
  for(i=0;i<A->r;i++)
   {
    for(j=0;j<A->c;j++)
     printf("%lu ",nmod_mat_entry(A,i,j));
    printf("\n");
   }
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
    // Should I use a GMP subroutine instead n_xgcd?
    g=n_xgcd(v,u,y,x);
    #if BUG0_nmod_mat_HNF
     printf("DKryskov_gcd_ui() po 0: u=%lu v=%lu\n",*u,*v);
    #endif
    assert( g == (*v)*y - (*u)*x );
    *u = n-( (*u) % n );
    *v %= n;
    return g;
   }
  g=n_xgcd(u,v,x,y);
    #if BUG0_nmod_mat_HNF
     printf("DKryskov_gcd_ui() po 1: u=%lu v=%lu\n",*u,*v);
    #endif
  assert( g == (*u)*x - (*v)*y );
  *v = n-( (*v) % n );
  *u %= n;
  return g;
 }

static __inline__ mp_limb_t
DKryskov_gcd_ui_2arg( mp_limb_t alpha, mp_limb_t betta )
 {
  if(alpha>betta)
   return n_gcd( alpha , betta );
  return n_gcd( betta , alpha );
 }

#define MPLUS( a, b, mod, n ) _nmod_add(a,b,mod)
#define MMUL( a, b, mod, n ) nmod_mul(a,b,mod)

static __inline__ void
DKryskov_nmod_zero_line(nmod_mat_t A,long i,long j,mp_limb_t* scrth)
 {
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_zero_line() i=%ld j=%ld\n",i,j);
  #endif
  assert(i != j);
  mp_limb_t x=nmod_mat_entry(A,i,i);
  mp_limb_t y=nmod_mat_entry(A,j,i);
  assert(x);
  assert(y);
  mp_limb_t n=A->mod.n;
  mp_limb_t u,v,g;
  assert(g);
  g=DKryskov_gcd_ui(&u,&v,x,y,n);
  #if 0
   assert( g == MPLUS( MMUL(u,x,A->mod,n), MMUL(v,y,A->mod,n), A->mod,n) );
   // somehow operands are correct, but sum is wrong, WTF? compiler error? 
   //MPLUS and MMUL only used in assert
  #endif
  #if 0 && !defined(NDEBUG)
   if( g != MPLUS( MMUL(u,x,A->mod,n), MMUL(v,y,A->mod,n), A->mod,n) )
    {
     x=17, y=7597145953435345532, u=9831600645622211865, v=12,
     u*x mod n = 1 - correct
     v*y mod n = 0 - correct
     their sum should be one, but expression MPLUS(...) evaluates to 
      3252452166838860553 instead, assert fails. WTF?
     printf("x=%lu, y=%lu\n",x,y);
     printf("n=%lu\n",n);
     printf("g=%lu\n",g);
     printf("u*x+v*y mod n=%lu\n", 
      MPLUS( MMUL(u,x,A->mod,n), MMUL(v,y,A->mod,n), A->mod,n) );
     printf("u=%lu, v=%lu\n",u,v);
     
     printf( "u*x mod n = %lu\n",MMUL(u,x,A->mod,n) );
     printf( "v*y mod n = %lu\n",MMUL(v,y,A->mod,n) );
     printf( "their sum again %lu\n", 
      MPLUS( MMUL(u,x,A->mod,n), MMUL(v,y,A->mod,n), A->mod,n) );
     
     assert(0);
    }
  #endif
  mp_limb_t iPLUS=i+1;
  mp_limb_t* alpha=A->rows[i]+iPLUS;
  mp_limb_t* betta=A->rows[j]+iPLUS;
  mp_limb_t vec_len=A->c-iPLUS;
  _nmod_vec_scalar_mul_nmod(    scrth, alpha, vec_len, u  , A->mod );
  _nmod_vec_scalar_addmul_nmod( scrth, betta, vec_len, v  , A->mod );
  _nmod_vec_scalar_mul_nmod(    betta, betta, vec_len, x/g , A->mod );
  _nmod_vec_scalar_addmul_nmod( betta, alpha, vec_len, n-y/g, A->mod );
  memcpy(alpha,scrth,vec_len*sizeof(mp_limb_t));
  alpha[-1]=g;
  betta[-1]=0;
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_zero_line() done\n");
   nmod_mat_print(A);
  #endif
 }

static __inline__ void 
DKryskov_nmod_Gauss_upper_last_col(nmod_mat_t A,long last_col)
 {
  mp_limb_t s=nmod_mat_entry( A, last_col, last_col );
  assert(s);
  long j;
  for(j=last_col;j--;)
   nmod_mat_entry( A, j, last_col ) %= s;
 }

static __inline__ void 
DKryskov_nmod_Gauss_upper(nmod_mat_t A)
 {
  long last_i=A->c-1;
  long j,i;
  mp_limb_t n=A->mod.n;
  mp_limb_t* sP,*tP;
  mp_limb_t s,t;
  for(i=1;i<last_i;i++)
   { 
    sP=&nmod_mat_entry( A, i,i );
    mp_limb_t s=*sP;
    assert(s);
    long v_len=A->c-i;
    for(j=i;j--;)
     {
      tP=&nmod_mat_entry( A, j,i );
      t=*tP;
      if(t)
       {
        mp_limb_t q=t/s;
        if(q)
         _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-q, A->mod );
        assert( *tP < s);
       }
     }
   }
  DKryskov_nmod_Gauss_upper_last_col( A, last_i );
 }

static __inline__ void 
DKryskov_nmod_early_abort(nmod_mat_t A,long e)
//all columns after e have 1 on diagonal, so save some effort
 {
  long m=A->c;
  long i,j;
  mp_limb_t* tP,* sP;
  long ePLUS=e+1;
  j=sizeof(mp_limb_t)*(m-ePLUS);
  for(i=m;i--;)
   {
    tP=A->rows[i]+ePLUS;
    memset(tP,0,j);
    #if BUG0_nmod_mat_HNF
     if(i>=ePLUS)
      {
       printf("putting 1, i=%ld, ePLUS=%ld\n",i,ePLUS);
      }
    #endif
    if(i>=ePLUS)
     tP[i-ePLUS]=1;
   }
    #if BUG0_nmod_mat_HNF
     printf("nmod_early_abort(), right half nulled\n");
     nmod_mat_print(A);
     printf("\n");
    #endif
  // now like Gauss_upper, but don't change anything after column e
  mp_limb_t n=A->mod.n;
  for(i=1;i<e;i++)
   { 
    sP=&nmod_mat_entry( A, i, i );
    mp_limb_t s=*sP;
    assert(s);
    long v_len=ePLUS-i;
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
    #if BUG0_nmod_mat_HNF
     printf("nmod_early_abort(), all but col no. %ld\n",e);
     nmod_mat_print(A);
     printf("\n");
    #endif
  DKryskov_nmod_Gauss_upper_last_col( A, e );
    #if BUG0_nmod_mat_HNF
     printf("nmod_early_abort() done\n",e);
     nmod_mat_print(A);
     printf("\n");
    #endif
 }

static __inline__ void
DKryskov_nmod_reduce_diag(nmod_mat_t A,long i,mp_limb_t det_tgt,mp_limb_t* t)
 {
  assert(i<A->c-1);
  if( det_tgt % nmod_mat_entry(A,i,i) )
   {
          #if BUG0_nmod_mat_HNF
           printf("nmod_mat_HNF(): gonna fix diagonal, i=%ld, new det_tgt=%ld\n",i,det_tgt);
           nmod_mat_print(A);
           printf("\n");
          #endif
    assert( 0 == nmod_mat_entry(A,i+1,i) );
    // Read DomichKannanTrotter87.pdf before asking me questions
    nmod_mat_entry(A,i+1,i)=det_tgt;
    DKryskov_nmod_zero_line(A,i,i+1,t);
   }
 }

static __inline__ void
DKryskov_nmod_reduce_last( mp_limb_t* se_corner, mp_limb_t det_tgt )
 {
    #if BUG0_nmod_mat_HNF
     printf("nmod_mat_HNF(): bad last element %ld, det_tgt=%ld\n",*se_corner,
      det_tgt);
    #endif
  *se_corner = DKryskov_gcd_ui_2arg( *se_corner, det_tgt );
 }

void
nmod_mat_HNF(nmod_mat_t A)
/*
HNF: Hermite Normal Form as defined by W.Stein/C.Pernet

Input data:
~~~~~~~~~~
A: square, non-singular over integers, A.r>0, det(A) is a multiple of n=A->mod.n
   (can be equal)

Output data:
~~~~~~~~~~~
A modified so new A = HNF of old A over Z
 
Literature:
~~~~~~~~~~
"Hermite Normal Form Computation Using Modulo Determinant Arithmetic"
aka DomichKannanTrotter87.pdf 
*/
 {
  long m=A->c;
  long i;
  mp_limb_t* scratch=(mp_limb_t*)malloc( sizeof(mp_limb_t) * (m-1) );
  long det_tgt=A->mod.n;
  for(i=0;i<m;i++)
   {
    // main loop: zap lower part, fix diagonal
    #if BUG0_nmod_mat_HNF
     printf("nmod_mat_HNF() i=%ld\n",i);
     nmod_mat_print(A);
     printf("\n");
    #endif
    long j=DKryskov_nmod_find_nonzero(A,i,det_tgt);
    #if BUG0_nmod_mat_HNF
     printf("nmod_mat_HNF() find_nonzero =%ld\n",j);
     printf("\n");
    #endif
    if(j==-2)
     {
      nmod_mat_entry(A,i,i)=det_tgt;
      if(i<m-1)
       {
        free(scratch);
        DKryskov_nmod_early_abort(A,i);
        return;
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
          // bad input: modulo not a multiple of matrice determinant
          DKryskov_nmod_reduce_last( &nmod_mat_entry(A,i,i), det_tgt );
          break;
         }
       }
      if( i==m-1 )
       //skip to Gauss_upper
       break; 
     }
    if(j>=0)
     {
      // some rows operations required
      DKryskov_nmod_zero_line(A,i,j,scratch);
      j=i;
      while(1)
       {
        j += 1;
        if(j>=m)
         break;
        if( nmod_mat_entry(A,j,i) )
         {
          DKryskov_nmod_zero_line(A,i,j,scratch);
          #if BUG0_nmod_mat_HNF
           printf("nmod_mat_HNF(): zapped line j=%ld with i=%ld\n",j,i);
           nmod_mat_print(A);
           printf("\n");
          #endif
         }
       }
     }
    DKryskov_nmod_reduce_diag(A,i,det_tgt,scratch);
    det_tgt /= nmod_mat_entry(A,i,i);
          #if BUG0_nmod_mat_HNF
           printf("nmod_mat_HNF(): fixed diagonal, i=%ld, new det_tgt=%ld\n",i,det_tgt);
           nmod_mat_print(A);
           printf("\n");
          #endif
    if( (det_tgt==1) && (i<m-1) )
     {
      free(scratch);
      DKryskov_nmod_early_abort(A,i);
      return;
     }
   }// main loop end
  free(scratch);
  DKryskov_nmod_Gauss_upper(A);
 }

#undef MPLUS
#undef MMUL
