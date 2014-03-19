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
   
  mod_mat_HNF.c should be easier to understand than this file
#endif

static __inline__ long 
DKryskov_nmod_find_nonzero(nmod_mat_t A,long col,mp_limb_t det_tgt)
/*
 -3 if diagonal element is 1
 -2 if all is zero
 -1 if single non-zero
 line>0 if A[line,i] is non-zero
*/
 {
  long i0=-1,i1=-1,k=col;
  mp_limb_t m0=det_tgt;
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
  printf("\n");
 }

static __inline__ void 
vec_print( const char* m, mp_limb_t* v, long s )
 {
  printf(m);
  long i;
  for(i=0;i<s;i++)
   printf("%lu ",v[i]);
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

static __inline__ long
DKryskov_nmod_zero_line(nmod_mat_t A,long i,long j,mp_limb_t n,mp_limb_t* scrth)
//returns 1 iff new A[i,i] becomes 1
 {
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_zero_line() i=%ld j=%ld\n",i,j);
  #endif
  assert(i != j);
  long iPLUS=i+1;
  mp_limb_t* alpha=A->rows[i]+iPLUS;
  mp_limb_t* betta=A->rows[j]+iPLUS;
  mp_limb_t x=alpha[-1];
  mp_limb_t y=betta[-1];
  assert(x>1);
  assert(y>1);
  mp_limb_t u,v,g;
  g=DKryskov_gcd_ui(&u,&v,x,y,n);
  assert(g);
  #if 0
   //MPLUS and MMUL only used in assert, assert below was seen failing
   assert( g == MPLUS( MMUL(u,x,A->mod,n), MMUL(v,y,A->mod,n), A->mod,n) );
   // somehow operands are correct, but sum is wrong, WTF? compiler error? 
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
  mp_limb_t vec_len=A->c-iPLUS;
  _nmod_vec_scalar_mul_nmod(    scrth, alpha, vec_len, u  , A->mod );
  _nmod_vec_scalar_addmul_nmod( scrth, betta, vec_len, v   , A->mod );
  _nmod_vec_scalar_mul_nmod(    betta, betta, vec_len, (x/g) , A->mod );
  _nmod_vec_scalar_addmul_nmod( betta, alpha, vec_len, (n-y/g) , A->mod );
  memcpy(alpha,scrth,vec_len*sizeof(mp_limb_t));
  alpha[-1]=g;
  betta[-1]=0;
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_zero_line() done\n");
   nmod_mat_print(A);
  #endif
  return g==1;
 }

static __inline__ void
DKryskov_nmod_easy_zl(nmod_mat_t A,long i,long j,mp_limb_t n)
//act like DKryskov_nmod_zero_line, but return nothing and use the fact that
// betta=0
 {
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_easy_zl() starts, i=%ld j=%ld n=%ld\n",i,j,n);
   nmod_mat_print(A);
  #endif
  assert(i != j);
  long iPLUS=i+1;
  mp_limb_t* alpha=A->rows[i]+iPLUS;
  mp_limb_t* betta=A->rows[j]+iPLUS;
  mp_limb_t x=alpha[-1];
  mp_limb_t y=betta[-1];
  assert(x>1);
  assert(y>1);
  mp_limb_t u,v,g;
  g=DKryskov_gcd_ui(&u,&v,x,y,n);
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_easy_zl(): x=%ld y=%ld u=%ld v=%ld g=%ld\n",x,y,u,v,g);
  #endif
  assert(g);
  mp_limb_t vec_len=A->c-iPLUS;
  _nmod_vec_scalar_mul_nmod( betta, alpha, vec_len, (n-y/g) , A->mod );
  _nmod_vec_scalar_mul_nmod( alpha, alpha, vec_len, u  , A->mod );
  alpha[-1]=g;
  betta[-1]=0;
  #if BUG0_nmod_mat_HNF
   printf("DKryskov_nmod_easy_zl() ends\n");
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
/*
zap A over diagonal
attempt to avoid row operations if possible
*/
 {
  long last_i=A->c-1;
  long j,i;
  mp_limb_t n=nmod_mat_entry(A,1,1);
  for(i=2;i<=last_i;i++)
   n *= nmod_mat_entry(A,i,i);
  #if BUG0_nmod_mat_HNF
   printf( "Gauss_upper() start, n=%lu\n",n );
  #endif
  mp_limb_t* sP,*tP;
  mp_limb_t s,t_ori,t_upd;
  for(i=1;i<last_i;i++)
   {
    sP=&nmod_mat_entry( A, i,i );
    mp_limb_t s=*sP;
    assert(s);
    long v_len=A->c-i;
    for(j=i;j--;)
     {
      tP=&nmod_mat_entry( A, j,i );
      t_ori=*tP;
      if(t_ori >= s)
       {
        t_upd=t_ori % n;
        if(t_upd < s)
         // don't need vector operation
         nmod_mat_entry( A, j,i ) = t_upd;
        else
         {
          mp_limb_t q=t_upd/s;
          #if BUG0_nmod_mat_HNF
           printf("subtracting row no. %ld from %ld, q=%lu\n",i,j,q);
          #endif
          _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-q, A->mod );
          *tP %= n;
          #if BUG0_nmod_mat_HNF
           nmod_mat_print(A);
          #endif
         }
       }
      assert( *tP < s );
     }
    n /= s; 
    #if BUG0_nmod_mat_HNF
     printf( "Gauss_upper() done with col %ld, new n=%lu\n" ,i ,n );
     nmod_mat_print(A);
    #endif
   }
  DKryskov_nmod_Gauss_upper_last_col( A, last_i );
 }

static __inline__ void 
DKryskov_nmod_early_abort(nmod_mat_t A,long e)
//all columns after e have 1 on diagonal, so save some effort
// TODO: apply decreasing modulo to avoid row operations
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
    #endif
  DKryskov_nmod_Gauss_upper_last_col( A, e );
    #if BUG0_nmod_mat_HNF
     printf("nmod_early_abort() done\n",e);
     nmod_mat_print(A);
    #endif
 }

static __inline__ void
DKryskov_nmod_reduce_diag(nmod_mat_t A,long i,mp_limb_t det_tgt,mp_limb_t* scratch)
 {
  assert(i<A->c-1);
  if( det_tgt % nmod_mat_entry(A,i,i) )
   {
          #if BUG0_nmod_mat_HNF
           printf("nmod_mat_HNF(): gonna fix diagonal, i=%ld, new det_tgt=%ld\n",i,det_tgt);
           nmod_mat_print(A);
          #endif
    assert( 0 == nmod_mat_entry(A,i+1,i) % det_tgt );
    // Read DomichKannanTrotter87.pdf before asking me questions
    nmod_mat_entry(A,i+1,i)=det_tgt;
    (void)DKryskov_nmod_zero_line(A,i,i+1,det_tgt,scratch);
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

static __inline__ void
DKryskov_nmod_1_lower(nmod_mat_t A,long col,long j,mp_limb_t n)
/*
A[col,col] is known to be 1
zap column col starting from row j
operate modulo n
*/
 {
    #if BUG0_nmod_mat_HNF
     printf("DKryskov_nmod_1_lower() start col=%ld j=%ld n=%lu\n",col,j,n);
     nmod_mat_print(A);
    #endif
  long m=A->c;
  mp_limb_t* sP=A->rows[col]+col;
  assert( 1 == *sP );
  long v_len = m-col;
  while(j < m)
   {
    mp_limb_t* tP=A->rows[j]+col;
    mp_limb_t t_ori=*tP;
    if(t_ori)
     {
      mp_limb_t t_upd=t_ori % n;
      if(t_upd)
       {
        #if BUG0_nmod_mat_HNF
         printf("q=%lu old vector=",n-t_upd);
         vec_print( "", tP, v_len );
        #endif
        _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-t_upd, A->mod );
        #if BUG0_nmod_mat_HNF
         vec_print( "new vector=", tP, v_len );
         printf("\n");
        #endif
       }
      else
       *tP = 0;
     }
    ++j;
   }
    #if BUG0_nmod_mat_HNF
     printf("DKryskov_nmod_1_lower() end\n");
     nmod_mat_print(A);
    #endif
 }

void
nmod_mat_HNF(nmod_mat_t A)
/*
HNF: Hermite Normal Form as defined by W.Stein/C.Pernet

Input data:
~~~~~~~~~~
A: square, non-singular over integers, n=A->mod.n is a multiple of abs(det(A))
  (can be equal)

Output data:
~~~~~~~~~~~
Define H := HNF( input A )

Output A contains matrice H: its diagonal and above entries match those of H:

for i<=j A[i,j]=H[i,j]
 
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
  for(i=0;;i++)
   {// main loop: zap below diagonal, fix diagonal, maintain decreasing modulo
    assert( (i>=0) && (i<m) );
    #if BUG0_nmod_mat_HNF
     printf("nmod_mat_HNF() i=%ld n=%lu\n",i,det_tgt);
     nmod_mat_print(A);
    #endif
    long j=DKryskov_nmod_find_nonzero(A,i,det_tgt);
    #if BUG0_nmod_mat_HNF
     printf("nmod_mat_HNF() find_nonzero =%ld\n",j);
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
      while(1)
       {
        if( ++j >= m )
         break;
        if( nmod_mat_entry(A,j,i) )
         {
          if(DKryskov_nmod_zero_line(A,i,j,det_tgt,scratch))
           {
            if( ++j < m )
             DKryskov_nmod_1_lower(A,i,j,det_tgt);
            break;
           }
          #if BUG0_nmod_mat_HNF
           printf("nmod_mat_HNF(): zapped line j=%ld with i=%ld\n",j,i);
           nmod_mat_print(A);
          #endif
         }
       }
     }
    DKryskov_nmod_reduce_diag(A,i,det_tgt,scratch);
    assert( 0 == det_tgt % nmod_mat_entry(A,i,i) );
    det_tgt /= nmod_mat_entry(A,i,i);
          #if BUG0_nmod_mat_HNF
           printf("nmod_mat_HNF(): fixed diagonal, i=%ld, new det_tgt=%ld\n",i,det_tgt);
           nmod_mat_print(A);
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
