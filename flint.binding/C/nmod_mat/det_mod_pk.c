// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// 4x4 determinant calculation is based on comments in Sage file
//  matrix_integer_dense.pyx

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#if WANT_ASSERT_IN_DET_MOD_PK
 #include <assert.h>
 #define ASSERT(x) assert(x)
#else
 #define ASSERT(x)
#endif
#include "../ulong_extras/ulong_extras_.h"
#include "../nmod_mat/nmod_mat_.h"

#define VECTOR_DOT_IN_cutoff_4 0
#define LOUD_4x4_INVERT 0
#define BUG_IN_cutoff_4 0
#define DUMP_cutoff_1_CALL 0
#if LOUD_4x4_INVERT
 #define PRINTF flint_printf
 #define SHOW_0TH_COL(m,a) flint_printf("%s 0th col: %wu %wu %wu %wu %wu\n", \
  m,nmod_mat_entry(a,0,0),nmod_mat_entry(a,1,0),nmod_mat_entry(a,2,0),      \
  nmod_mat_entry(a,3,0),nmod_mat_entry(a,4,0));
#else
 #define PRINTF
 #define SHOW_0TH_COL
#endif

#if GMP_LIMB_BITS == 64 && defined (__amd64__) && 1
 #define ALIGN_inv_array 0x10
 #include <xmmintrin.h> // need __m128i
 #undef uint128_t
 #define uint128_t __m128i
#endif

#if defined(WX_MINUS_YZ) && 1
 #define WX_MINUS_YZ_5arg( rez, w,x,y,z ) \
  WX_MINUS_YZ( rez, w,x,y,z, n,i );
#else
 #define WX_MINUS_YZ_5arg( rez, w,x,y,z ) \
  rez=SUB_mod_n( MUL_mod_n(w,x), MUL_mod_n(y,z) );
#endif

#if defined(MULADD_pk) && 1
 #define MULADD_3arg(r,a,b) MULADD_pk(r, a,b, n,i);
#else
 #define MULADD_3arg(r,a,b) r=ADD_mod_n(r, MUL_mod_n(a,b) );
#endif

//TODO: use better algorithm for vector dot product, in all subroutines

#define MUL n_mulmod_preinv_4arg
#define MUL_mod_n(x,y) n_mulmod_preinv_4arg(x,y,n,i)
#define SUB_mod_n(x,y) n_submod(x,y,n)
#define ADD_mod_n(x,y) n_addmod(x,y,n)
#if defined(VECTOR_DOT_HEAD) && 0
 #if SPEEDUP_NMOD_RED3
 #else
  #error "can't compile this: need 2**128 modulo n"
 #endif
 #define DET_4x4                               \
  WX_MINUS_YZ_5arg(a, r0[3],r1[2], r0[2],r1[3]) \
  WX_MINUS_YZ_5arg(b, r2[1],r3[0], r2[0],r3[1]) \
  VECTOR_DOT_HEAD(a,b);                         \
  WX_MINUS_YZ_5arg(a, r0[1],r1[3], r0[3],r1[1]) \
  WX_MINUS_YZ_5arg(b, r2[2],r3[0], r2[0],r3[2]) \
  VECTOR_DOT_BODY(a,b);                         \
  WX_MINUS_YZ_5arg(a, r0[2],r1[1], r0[1],r1[2]) \
  WX_MINUS_YZ_5arg(b, r2[3],r3[0], r2[0],r3[3]) \
  VECTOR_DOT_BODY(a,b);                         \
  WX_MINUS_YZ_5arg(a, r0[3],r1[0], r0[0],r1[3]) \
  WX_MINUS_YZ_5arg(b, r2[2],r3[1], r2[1],r3[2]) \
  VECTOR_DOT_BODY(a,b);                         \
  WX_MINUS_YZ_5arg(a, r0[0],r1[2], r0[2],r1[0]) \
  WX_MINUS_YZ_5arg(b, r2[3],r3[1], r2[1],r3[3]) \
  VECTOR_DOT_BODY(a,b);                         \
  WX_MINUS_YZ_5arg(a, r0[1],r1[0], r0[0],r1[1]) \
  WX_MINUS_YZ_5arg(b, r2[3],r3[2], r2[2],r3[3]) \
  VECTOR_DOT_BODY(a,b);                         \
  VECTOR_DOT_TAIL(r, n,i,M->mod.norm )
#else
 #define DET_4x4                               \
  WX_MINUS_YZ_5arg(a, r0[3],r1[2], r0[2],r1[3]) \
  WX_MINUS_YZ_5arg(b, r2[1],r3[0], r2[0],r3[1]) \
  r=MUL_mod_n( a, b );                          \
  WX_MINUS_YZ_5arg(a, r0[1],r1[3], r0[3],r1[1]) \
  WX_MINUS_YZ_5arg(b, r2[2],r3[0], r2[0],r3[2]) \
  MULADD_3arg(r, a,b)                           \
  WX_MINUS_YZ_5arg(a, r0[2],r1[1], r0[1],r1[2]) \
  WX_MINUS_YZ_5arg(b, r2[3],r3[0], r2[0],r3[3]) \
  MULADD_3arg(r, a,b)                           \
  WX_MINUS_YZ_5arg(a, r0[3],r1[0], r0[0],r1[3]) \
  WX_MINUS_YZ_5arg(b, r2[2],r3[1], r2[1],r3[2]) \
  MULADD_3arg(r, a,b)                           \
  WX_MINUS_YZ_5arg(a, r0[0],r1[2], r0[2],r1[0]) \
  WX_MINUS_YZ_5arg(b, r2[3],r3[1], r2[1],r3[3]) \
  MULADD_3arg(r, a,b)                           \
  WX_MINUS_YZ_5arg(a, r0[1],r1[0], r0[0],r1[1]) \
  WX_MINUS_YZ_5arg(b, r2[3],r3[2], r2[2],r3[3]) \
  MULADD_3arg(r, a,b)
#endif

static __inline__ mp_limb_t
nmod_mat_det_dim4(const nmod_mat_t M)
/*
count determinant of 4x4 submatrice in upper-left corner
 
Use same algorithm as Sage:
 a = x[3]*x[6]; a-=x[2]*x[7]; b = x[ 9]*x[12]; b-=x[8]*x[13];  r  = a*b
 a = x[1]*x[7]; a-=x[3]*x[5]; b = x[10]*x[12]; b-=x[8]*x[14];  r += a*b
 a = x[2]*x[5]; a-=x[1]*x[6]; b = x[11]*x[12]; b-=x[8]*x[15];  r += a*b
 a = x[3]*x[4]; a-=x[0]*x[7]; b = x[10]*x[13]; b-=x[9]*x[14];  r += a*b
 a = x[0]*x[6]; a-=x[2]*x[4]; b = x[11]*x[13]; b-=x[9]*x[15];  r += a*b
 a = x[1]*x[4]; a-=x[0]*x[5]; b = x[11]*x[14]; b-=x[10]*x[15]; r += a*b
*/
 {
  // initialization from incompatible pointer type for next line --- why?
  const mp_limb_t** const rows=M->rows;
  const mp_limb_t* r0=rows[0];
  const mp_limb_t* r1=rows[1];
  const mp_limb_t* r2=rows[2];
  const mp_limb_t* r3=rows[3];
  const mp_limb_t n=M->mod.n;
  const mp_limb_t i=M->mod.ninv;
  mp_limb_t r,a,b;
  DET_4x4
  return r;
 }

static __inline__ mp_limb_t
nmod_mat_det_dim4_SE(const nmod_mat_t M)
/*
count determinant of 4x4 submatrice in lower-right corner
use same algorithm as nmod_mat_det_dim4()
*/
 {
  slong s=M->r-4;
  const mp_limb_t** const rows=M->rows;
  const mp_limb_t* r0=rows[s]+s;
  const mp_limb_t* r1=rows[s+1]+s;
  const mp_limb_t* r2=rows[s+2]+s;
  const mp_limb_t* r3=rows[s+3]+s;
  const mp_limb_t n=M->mod.n;
  const mp_limb_t i=M->mod.ninv;
  mp_limb_t r,a,b;
  DET_4x4
  return r;
 }

#if defined(WX_MINUS_YZ)
 #define DIM2_DET( rez, r0, r1 ) \
  WX_MINUS_YZ( rez, r0[0],r1[1], r0[1],r1[0], n,i );
#else
 #define DIM2_DET( rez, r0, r1 )          \
   rez=          MUL( r0[0], r1[1], n, i ); \
   rez=n_submod( rez,                       \
                 MUL( r0[1], r1[0], n, i ), \
             n );
#endif

static __inline__ mp_limb_t
nmod_mat_det_dim3(const nmod_mat_t A)
 {
  const mp_limb_t** const rows=A->rows;
  const mp_limb_t* r0=rows[0];
  const mp_limb_t* r1=rows[1];
  const mp_limb_t* r2=rows[2];
  const mp_limb_t n=A->mod.n;
  const mp_limb_t i=A->mod.ninv;
  mp_limb_t rez,t;
  #if defined(VECTOR_DOT_HEAD) && defined(SPEEDUP_NMOD_RED3)
   DIM2_DET( t, r0, r1 );
   VECTOR_DOT_HEAD( t, r2[2] );
   DIM2_DET( t, r1, r2 );
   VECTOR_DOT_BODY( t, r0[2] );
   rez=r1[2];
   if(rez)
    {
     rez=n-rez;
     DIM2_DET( t, r0, r2 );
     VECTOR_DOT_BODY( t, rez );
    }
   VECTOR_DOT_TAIL( t, n,i,A->mod.norm );
   return t;
  #else
   DIM2_DET( rez, r0, r1 );
   rez=MUL( rez, r2[2], n, i );
   DIM2_DET( t, r0, r2 );
   t=MUL( t, r1[2], n, i );
   rez=n_submod( rez, t, n);
   DIM2_DET( t, r1, r2 );
   t=MUL( t, r0[2], n, i );
   return n_addmod( rez, t, n);
  #endif
 }

static __inline__ mp_limb_t
nmod_mat_det_dim2(const nmod_mat_t A)
 {
  const mp_limb_t* r0=A->rows[0];
  const mp_limb_t* r1=A->rows[1];
  const mp_limb_t n=A->mod.n;
  const mp_limb_t i=A->mod.ninv;
  mp_limb_t rez;
  DIM2_DET( rez, r0, r1 );
  return rez;
 }

#undef DIM2_DET
#undef MUL
#undef MUL_mod_n
#undef ADD_mod_n
#undef SUB_mod_n
#undef DET_4x4
#undef WX_MINUS_YZ_5arg
#undef MULADD_3arg

static __inline__ mp_limb_t
det_mod_pk_SE_0th_row(nmod_mat_t M,slong* negate_det,mp_limb_t p)
/*
arrange it so M[-1,-1] is non-degenerate
return the pivot element on success, 0 on failure
*/
 {
  slong dim=M->r;
  slong s=dim-1,i;
  mp_limb_t found;
  mp_limb_t** rows=M->rows;
  if( (found=rows[s][s]) % p )
   return found;
  for(i=s;i--;)
   {
    if( (found=rows[i][s]) % p )
     {
      (*negate_det) ^= 1;
      MP_PTR_SWAP( rows[s], rows[i] );
      return found;
     }
   }
  return 0;
 }

static __inline__ void 
det_mod_pk_SE_2x2_invert(mp_limb_t* invM,nmod_mat_t M,mp_limb_t alpha_inv,mp_limb_t betta,
  mp_limb_t epsln_inv, mp_limb_t gamma)
/*
X=2x2 matrice in lower-right corner of M,

    ( delta gamma )   ( 1  gamma*alpha' )   ( epsln   0    )
X = (             ) = (                 ) * (              )
    ( betta alpha )   ( 0       1       )   ( betta  alpha )
 
Invert X using formula

 '   ( epsln'   0    )   ( 1 -gamma*alpha' )
X  = (               ) * (                 )
     ( zetta  alpha' )   ( 0       1       )

where zetta=-epsln'*alpha'*betta
 
Store inverted X into SE corner of invM. 
Produce result modulo M->mod.n.
*/
 {
  const nmod_t mod=M->mod;
  mp_limb_t zetta=mod.n - epsln_inv;
  zetta=n_mulmod_preinv_4arg(zetta, alpha_inv, mod.n, mod.ninv);
  invM[14]=zetta=n_mulmod_preinv_4arg(zetta, betta    , mod.n, mod.ninv);
  mp_limb_t tempp=mod.n - alpha_inv;
  tempp=n_mulmod_preinv_4arg(tempp, gamma, mod.n, mod.ninv);
  invM[10]=epsln_inv;
  invM[11]=n_mulmod_preinv_4arg( epsln_inv, tempp, mod.n, mod.ninv);
  //invM[14]=zetta;
  #if 0
   tempp=n_mulmod_preinv_4arg(tempp, zetta, mod.n, mod.ninv);
   invM[15]=n_addmod( tempp, alpha_inv, mod.n );
  #else
   MULADD_pk( alpha_inv, tempp,zetta, mod.n, mod.ninv );
   invM[15]=alpha_inv;
  #endif
  #if LOUD_4x4_INVERT
   flint_printf("SE 2x2 corner inverse: %wu %wu | %wu %wu\n",
    invM[10],invM[11],invM[14],invM[15]);
  #endif
 }

static __inline__ mp_limb_t 
det_mod_pk_SE_1st_row(mp_limb_t* invM,nmod_mat_t M,slong* negate_det,
  const p_k_pk_t pp)
// return 1 on success  
 {
  mp_limb_t** rows=M->rows;
  const mp_limb_t dim_minus_1=M->r-1;
  const mp_limb_t dim_minus_2=dim_minus_1-1;
  const nmod_t mod=M->mod;
  const mp_limb_t alpha=invM[15];
  const mp_limb_t alpha_inv=inv_mod_pk_3arg(alpha,pp,mod);
  invM[0]=alpha_inv;
  const mp_limb_t betta=rows[dim_minus_1][dim_minus_2];
  const mp_limb_t alpha_inv_by_beta=n_mulmod_preinv_4arg(alpha_inv,betta,
   mod.n, mod.ninv );
  mp_limb_t* tail;
  mp_limb_t gamma,delta,epsln;
  //flint_printf("p,alpha,alpha_inv=%wu,%wu,%wu\n",p,alpha,alpha_inv);
  slong i;
  for(i=dim_minus_1;i--;)
   {
    tail=rows[i]+dim_minus_2;
    gamma=tail[1];
    delta=tail[0];
    //flint_printf("delta,gamma=%wu,%wu\n",delta%p_deg_k,gamma%p_deg_k);
    epsln=n_mulmod_preinv_4arg( gamma, alpha_inv_by_beta, mod.n, mod.ninv );
    epsln=n_submod( delta, epsln, mod.n );
    //flint_printf("i=%ld epsilon=%wu\n",i,epsln%p_deg_k);
    if( delta=(epsln % pp.p) )
     {
      if(i != dim_minus_2)
       {
        (*negate_det) ^= 1;
        MP_PTR_SWAP( rows[i], rows[dim_minus_2] );
       }
      invM[12]=epsln; // put the pivot into SW corner of invM
      epsln=inv_mod_pk_4arg(epsln,delta,pp,mod);
      det_mod_pk_SE_2x2_invert( invM, M, alpha_inv, betta, epsln, gamma );
      return 1;
     }
   }
  return 0;
 }

#define MUL(x,y) n_mulmod_preinv_4arg(x,y,mod.n,mod.ninv)
#define ADD(x,y) n_addmod(x,y,mod.n)
#define SUB(x,y) n_submod(x,y,mod.n)
#if defined(MULADD_pk)
 #define MULADD_3arg(r,a,b) MULADD_pk(r, a,b, mod.n,mod.ninv)
#endif

static __inline__ void 
det_mod_pk_mul_2x2( 
  mp_limb_t* r0,mp_limb_t* r1,
  const mp_limb_t* a0,const mp_limb_t* a1,
  const mp_limb_t* b0,const mp_limb_t* b1,
  const nmod_t mod)
 {
  #if defined(VECTOR_DOT_2)
   VECTOR_DOT_2( r0[0], a0[0],b0[0], a0[1],b1[0], mod)
   VECTOR_DOT_2( r0[1], a0[0],b0[1], a0[1],b1[1], mod)
   VECTOR_DOT_2( r1[0], a1[0],b0[0], a1[1],b1[0], mod)
   VECTOR_DOT_2( r1[1], a1[0],b0[1], a1[1],b1[1], mod)
  #else
   #if defined(MULADD_3arg)
    mp_limb_t t;
    t=MUL(a0[0],b0[0]); MULADD_3arg(t, a0[1],b1[0]); r0[0]=t;
    t=MUL(a0[0],b0[1]); MULADD_3arg(t, a0[1],b1[1]); r0[1]=t;
    t=MUL(a1[0],b0[0]); MULADD_3arg(t, a1[1],b1[0]); r1[0]=t;
    t=MUL(a1[0],b0[1]); MULADD_3arg(t, a1[1],b1[1]); r1[1]=t;
   #else
    r0[0]=ADD( MUL(a0[0],b0[0]), MUL(a0[1],b1[0]) );
    r0[1]=ADD( MUL(a0[0],b0[1]), MUL(a0[1],b1[1]) );
    r1[0]=ADD( MUL(a1[0],b0[0]), MUL(a1[1],b1[0]) );
    r1[1]=ADD( MUL(a1[0],b0[1]), MUL(a1[1],b1[1]) );
   #endif
  #endif
 }

#undef MULADD_3arg

static __inline__ void 
det_mod_pk_mul_sub_2x2(
  mp_limb_t* eps0,mp_limb_t* eps1,   
  mp_limb_t* sou0,mp_limb_t* sou1,   
  mp_limb_t* zet0,mp_limb_t* zet1,
  const nmod_t mod)
/*
set epsilon=delta-gamma*zeta
delta is at sou0,sou1, gamma at sou0+2,sou1+2
*/
 {
  #if defined(VECTOR_DOT_2)
   // for non-inline version, 170 instructions
   mp_limb_t t;
   VECTOR_DOT_2(t, sou0[2],zet0[0], sou0[3],zet1[0], mod);
   eps0[0]=SUB( sou0[0], t );
   VECTOR_DOT_2(t, sou0[2],zet0[1], sou0[3],zet1[1], mod);
   eps0[1]=SUB( sou0[1], t );
   VECTOR_DOT_2(t, sou1[2],zet0[0], sou1[3],zet1[0], mod);
   eps1[0]=SUB( sou1[0], t );
   VECTOR_DOT_2(t, sou1[2],zet0[1], sou1[3],zet1[1], mod);
   eps1[1]=SUB( sou1[1], t );
  #else
   // for non-inline version, 142+8*22=318 instructions
   eps0[0]=SUB( sou0[0], ADD( MUL(sou0[2],zet0[0]), MUL(sou0[3],zet1[0]) ) );
   eps0[1]=SUB( sou0[1], ADD( MUL(sou0[2],zet0[1]), MUL(sou0[3],zet1[1]) ) );
   eps1[0]=SUB( sou1[0], ADD( MUL(sou1[2],zet0[0]), MUL(sou1[3],zet1[0]) ) );
   eps1[1]=SUB( sou1[1], ADD( MUL(sou1[2],zet0[1]), MUL(sou1[3],zet1[1]) ) );
  #endif
 }

#define ROW_23( i, j )                                          \
 cand_0=rows[i]+dim_minus_4;                                          \
 cand_1=rows[j]+dim_minus_4;                                              \
 det_mod_pk_mul_sub_2x2(invM,epsi1,   cand_0,cand_1,   zeta0,zeta1,   mod); \
 PRINTF("i=%w j=%w epsi=%wu %wu | %wu %wu\n",i,j,       \
  invM[0],invM[1],epsi1[0],epsi1[1]);                   \
 ok=SUB( MUL(invM[0],epsi1[1]), MUL(invM[1],epsi1[0]) );                     \
 if(0 == ok % p )                                                             \
  ok=0;

static __inline__ mp_limb_t 
det_mod_pk_SE_row_23(mp_limb_t* invM,nmod_mat_t M,slong* negate_det,mp_limb_t p)
/*
count 2x2 matrice epsilon=delta-gamma*alpha'*betta

                            (  delta gamma )
where SE corner of M equals (              )
                            (  betta alpha )

count det=determinant of epsilon, check if it is non-degenerate

return 0 on failure, det on success 
*/
 {
  const slong dim_minus_3 = M->r-3;
  const slong dim_minus_4 = dim_minus_3-1;
  mp_limb_t** rows=M->rows;
  const nmod_t mod=M->mod;
  // use NE corner to keep product alpha'*betta
  mp_limb_t* zeta0=invM+2;
  mp_limb_t* zeta1=invM+6;
  mp_limb_t* cand_0;
  mp_limb_t* cand_1;
  det_mod_pk_mul_2x2( zeta0,zeta1,   invM+10,invM+14,  
   rows[dim_minus_3+1]+dim_minus_4,rows[dim_minus_3+2]+dim_minus_4,   mod );
  #if LOUD_4x4_INVERT
   flint_printf("det_mod_pk_SE_row_23() zeta=%wu %wu | %wu %wu\n",
    zeta0[0],zeta0[1],zeta1[0],zeta1[1]);
  #endif
  //put epsilon into NW corner of invM
  //mp_limb_t* epsi0=invM+0; 
  mp_limb_t* epsi1=invM+4;
  slong i,j;
  mp_limb_t ok;
  ROW_23( dim_minus_4, dim_minus_3 );
  if( ok )
   return ok;
  ROW_23( 0, dim_minus_3 );
  if( ok )
   {
    (*negate_det) ^= 1;
    MP_PTR_SWAP( rows[0], rows[dim_minus_4] );
    return ok;
   }
  ROW_23( dim_minus_4, 0 );
  if( ok )
   {
    (*negate_det) ^= 1;
    MP_PTR_SWAP( rows[0], rows[dim_minus_3] );
    return ok;
   }
  //only try neighbors
  for(i=dim_minus_4-2;i>=0;i-=2)
   {
    j=i+1;
    ROW_23( i, j );
    if( ok )
     {
      MP_PTR_SWAP( rows[i], rows[dim_minus_4] );
      MP_PTR_SWAP( rows[j], rows[dim_minus_3] );
      return ok;
     }
   }
  return 0;
 }

#undef MUL
#undef ADD
#undef SUB
#undef ROW_23
#undef MULADD_3arg

static __inline__ mp_limb_t
det_mod_pk_examine_last_column(nmod_mat_t M,slong* negate_det,
  const p_k_pk_t pp)
/*
return 0 if last column is zero modulo p**k
otherwise find best element, put it into lower-right corner, return 1
*/
 {
  slong all_zero=1000;
  slong best_row=-1,i,j;
  const slong dim=M->r;
  const slong shi=dim-1;
  mp_limb_t** rows=M->rows;
  mp_limb_t e;
  // find less spoilt element
  for(i=dim;i--;)
   {
    // TODO: optimize for huge p
    e=rows[i][shi] % pp.p_deg_k;
    if(e)
     {
      j=n_remove( &e, pp.p );// TODO: maybe use n_remove2_precomp() here?
      if(j<all_zero)
       {
        all_zero=j;
        best_row=i;
        if( 1==j )
         break;
       }
     }
   }
  if(best_row == -1)
   return 0;
  if(best_row != shi)
   {
    (*negate_det) ^= 1;
    MP_PTR_SWAP( rows[best_row], rows[shi] );
   }
  return 1;
 }

static __inline__ mp_limb_t 
det_mod_pk_SE_corner_det(mp_limb_t* invM,nmod_mat_t M,mp_limb_t d)
//multiply d by alpha*epsilon where alpha is 0th pivot and epsilon is 1st pivot
 {
  const nmod_t mod=M->mod;
  const slong shi=M->r-1;
  PRINTF("det_mod_pk_SE_corner_det(): pivots=%wu %wu\n",M->rows[shi][shi],
   invM[12]);
  d=n_mulmod_preinv_4arg(d,M->rows[shi][shi],mod.n,mod.ninv);
  d=n_mulmod_preinv_4arg(d,invM[12]         ,mod.n,mod.ninv);
  return d;
 }

void
nmod_invert_2x2(mp_limb_t* s0,mp_limb_t* s1,const nmod_t mod,const p_k_pk_t pp)
 {
  mp_limb_t delta,gamma,betta,alpha=s1[1];
  mp_limb_t epsil,zetta,theta;
  // TODO: optimize for huge p
  if( betta=(alpha % pp.p) )
   {
    alpha=inv_mod_pk_4arg(alpha,betta,pp,mod);
    betta=s1[0];
    delta=s0[0];
    gamma=s0[1];
    // epsil=delta-gamma*alpha'*betta
    epsil=n_mulmod_preinv_4arg(gamma,alpha,mod.n,mod.ninv);
    epsil=n_mulmod_preinv_4arg(epsil,betta,mod.n,mod.ninv);
    epsil=n_submod(delta,epsil,mod.n);
    // zetta=-epsil'*alpha'*betta
    s0[0]=epsil=inv_mod_pk_3arg(epsil,pp,mod);
    ASSERT(epsil<pp.p_deg_k);
    zetta=n_mulmod_preinv_4arg(pp.p_deg_k-epsil,alpha,mod.n,mod.ninv);
    s1[0]=zetta=n_mulmod_preinv_4arg(zetta,betta,mod.n,mod.ninv);
    // theta=gamma*alpha'
    theta=n_mulmod_preinv_4arg(pp.p_deg_k-alpha,gamma,mod.n,mod.ninv);
    //r0[0]=epsil;
    s0[1]=n_mulmod_preinv_4arg(epsil,theta,mod.n,mod.ninv);
    //r1[0]=zetta;
    zetta=n_mulmod_preinv_4arg(zetta,theta,mod.n,mod.ninv);
    s1[1]=n_addmod(zetta,alpha,mod.n);
   }
  else
   {
    delta=s1[0];
    gamma=s0[1];
    betta=s0[0];
    // same as above, only switch gamma<->alpha and switch cols of result
    gamma=inv_mod_pk_3arg(gamma,pp,mod);
    // epsil=delta-gamma*alpha'*betta
    epsil=n_mulmod_preinv_4arg(gamma,alpha,mod.n,mod.ninv);
    epsil=n_mulmod_preinv_4arg(epsil,betta,mod.n,mod.ninv);
    epsil=n_submod(delta,epsil,mod.n);
    s0[1]=epsil=inv_mod_pk_3arg(epsil,pp,mod);
    ASSERT(epsil<pp.p_deg_k);
    zetta=n_mulmod_preinv_4arg(pp.p_deg_k-epsil,gamma,mod.n,mod.ninv);
    s1[1]=zetta=n_mulmod_preinv_4arg(zetta,betta,mod.n,mod.ninv);
    theta=n_mulmod_preinv_4arg(pp.p_deg_k-gamma,alpha,mod.n,mod.ninv);
    //r0[0]=epsil;
    s0[0]=n_mulmod_preinv_4arg(epsil,theta,mod.n,mod.ninv);
    //r1[0]=zetta;
    zetta=n_mulmod_preinv_4arg(zetta,theta,mod.n,mod.ninv);
    s1[0]=n_addmod(zetta,gamma,mod.n);
   }
 }

static __inline__ void
det_mod_pk_mul_negate_gamma_alphaINV( 
  mp_limb_t* gam0, mp_limb_t* gam1, 
  mp_limb_t* alf0, mp_limb_t* alf1, 
  const nmod_t mod)
 {
  mp_limb_t g00=n_negmod( gam0[0], mod.n );
  mp_limb_t g01=n_negmod( gam0[1], mod.n );
  mp_limb_t g10=n_negmod( gam1[0], mod.n );
  mp_limb_t g11=n_negmod( gam1[1], mod.n );
  #if defined(VECTOR_DOT_2)
   // non-inline det_mod_pk_mul_negate_gamma_alphaINV(): 170 instrucrions
   VECTOR_DOT_2(gam0[0], g00,alf0[0], g01,alf1[0], mod);
   VECTOR_DOT_2(gam0[1], g00,alf0[1], g01,alf1[1], mod);
   VECTOR_DOT_2(gam1[0], g10,alf0[0], g11,alf1[0], mod);
   VECTOR_DOT_2(gam1[1], g10,alf0[1], g11,alf1[1], mod);
  #else
   // non-inline det_mod_pk_mul_negate_gamma_alphaINV(): 121+8*22=297 instrucrions
   gam0[0] = n_addmod(
    n_mulmod_preinv_4arg( g00, alf0[0], mod.n, mod.ninv ),
    n_mulmod_preinv_4arg( g01, alf1[0], mod.n, mod.ninv ),
    mod.n );
   gam0[1] = n_addmod(
    n_mulmod_preinv_4arg( g00, alf0[1], mod.n, mod.ninv ),
    n_mulmod_preinv_4arg( g01, alf1[1], mod.n, mod.ninv ),
    mod.n );
   gam1[0] = n_addmod(
    n_mulmod_preinv_4arg( g10, alf0[0], mod.n, mod.ninv ),
    n_mulmod_preinv_4arg( g11, alf1[0], mod.n, mod.ninv ),
    mod.n );
   gam1[1] = n_addmod(
    n_mulmod_preinv_4arg( g10, alf0[1], mod.n, mod.ninv ),
    n_mulmod_preinv_4arg( g11, alf1[1], mod.n, mod.ninv ),
    mod.n );
  #endif
 }

void
det_mod_pk_mul_add_2x2(
  mp_limb_t* tgt0,mp_limb_t* tgt1,
  const mp_limb_t* rho0,const mp_limb_t* rho1,
  const mp_limb_t* eta0,const mp_limb_t* eta1,
  const nmod_t mod)
// tgt := tgt + rho * eta
 {
  #if defined(VECTOR_DOT_2_add)
   // 152 instructions
   register mp_limb_t t;
   #define DOT2(tgt, x0,y0, x1,y1) VECTOR_DOT_2_add(tgt, x0,y0, x1,y1, mod)
   DOT2(tgt0[0], rho0[0],eta0[0], rho0[1],eta1[0]);
   DOT2(tgt0[1], rho0[0],eta0[1], rho0[1],eta1[1]);
   DOT2(tgt1[0], rho1[0],eta0[0], rho1[1],eta1[0]);
   DOT2(tgt1[1], rho1[0],eta0[1], rho1[1],eta1[1]);
   #undef DOT2
  #else
   // 139+8*22=315 instructions
   tgt0[0] = n_addmod( tgt0[0],
     n_mulmod_preinv_4arg(rho0[0],eta0[0], mod.n,mod.ninv),
     mod.n
    );
   tgt0[0] = n_addmod( tgt0[0],
     n_mulmod_preinv_4arg(rho0[1],eta1[0], mod.n,mod.ninv),
     mod.n
    );
   tgt0[1] = n_addmod( tgt0[1],
     n_mulmod_preinv_4arg(rho0[0],eta0[1], mod.n,mod.ninv),
     mod.n
    );
   tgt0[1] = n_addmod( tgt0[1],
     n_mulmod_preinv_4arg(rho0[1],eta1[1], mod.n,mod.ninv),
     mod.n
    );
   tgt1[0] = n_addmod( tgt1[0],
     n_mulmod_preinv_4arg(rho1[0],eta0[0], mod.n,mod.ninv),
     mod.n
    );
   tgt1[0] = n_addmod( tgt1[0],
     n_mulmod_preinv_4arg(rho1[1],eta1[0], mod.n,mod.ninv),
     mod.n
    );
   tgt1[1] = n_addmod( tgt1[1],
     n_mulmod_preinv_4arg(rho1[0],eta0[1], mod.n,mod.ninv),
     mod.n
    );
   tgt1[1] = n_addmod( tgt1[1],
     n_mulmod_preinv_4arg(rho1[1],eta1[1], mod.n,mod.ninv),
     mod.n
    );
  #endif
 }

#define SHOW_VECTOR_4(i)                               \
 v=a->rows[i]+(a->r-4);                                 \
 flint_printf("%wu %wu %wu %wu\n",v[0]%n,v[1]%n,v[2]%n,v[3]%n);

#define SHOW_VEKTOR_4(v,n) \
 flint_printf("%wu %wu %wu %wu\n",(v)[0]%n,(v)[1]%n,(v)[2]%n,(v)[3]%n);

void
show_all_but_SE_4x4_corner(nmod_mat_t a,mp_limb_t n)
 {
  slong i,j;
  slong dim_minus_4=a->r-4;
  for(i=0;i<a->r;i++)
   {
    for(j=0;j<a->r;j++)
     {
      if( (i>=dim_minus_4) && (j>=dim_minus_4) )
       printf("* ");
      else
       flint_printf("%wu ",a->rows[i][j]);
     }
    printf("\n");
   }
 }

void 
show_SE_4x4_corner(nmod_mat_t a,mp_limb_t n)
 {
  mp_limb_t* v;
  SHOW_VECTOR_4(a->r-4)
  SHOW_VECTOR_4(a->r-3)
  SHOW_VECTOR_4(a->r-2)
  SHOW_VECTOR_4(a->r-1)
 }

void
show_Ainv(mp_limb_t* m,mp_limb_t n)
 {
  SHOW_VEKTOR_4(m,n);
  SHOW_VEKTOR_4(m+4,n);
  SHOW_VEKTOR_4(m+8,n);
  SHOW_VEKTOR_4(m+12,n);
 }

static __inline__ void
det_mod_pk_SE_4x4_invert(mp_limb_t* invM,nmod_mat_t M,const p_k_pk_t pp)
/*
Invert 4x4 matrice X stored at SE corner of M.

Put result into invM.

Use same formula as det_mod_pk_SE_2x2_invert(): 
for 
    ( delta gamma )   ( 1  gamma*alpha' )   ( epsln   0    )
X = (             ) = (                 ) * (              )
    ( betta alpha )   ( 0       1       )   ( betta  alpha )
 
inverse equals

 '   ( epsln'   0    )   ( 1 -gamma*alpha' )
X  = (               ) * (                 )
     ( zetta  alpha' )   ( 0       1       )

where zetta = -alpha'*betta*epsln', epsln=delta-gamma*alpha'*beta

alpha' is at SE corner of invM; epsln is at NW corner of invM
*/
 {
  #define p_deg_k pp.p_deg_k
  const nmod_t mod=M->mod;
  mp_limb_t** rows=M->rows;
  const slong dim_minus_2=M->r-2;
  #if LOUD_4x4_INVERT
   flint_printf("SE corner to be inverted:\n");
   show_SE_4x4_corner(M,p_deg_k);
  #endif
  // NW corner of target: epsilon->epsilon'
  PRINTF("det_mod_pk_SE_4x4_invert(): epsilon =%wu %wu | %wu %wu\n",
   invM[0]%p_deg_k,invM[1]%p_deg_k,invM[4]%p_deg_k,invM[5]%p_deg_k);
  nmod_invert_2x2(invM,invM+4,mod,pp);
  PRINTF("det_mod_pk_SE_4x4_invert(): epsilon'=%wu %wu | %wu %wu\n",
   invM[0]%p_deg_k,invM[1]%p_deg_k,invM[4]%p_deg_k,invM[5]%p_deg_k);
  // NE corner of source: gamma -> -gamma*alpha'
  mp_limb_t* source_NE_0=rows[dim_minus_2-2]+dim_minus_2;
  mp_limb_t* source_NE_1=rows[dim_minus_2-1]+dim_minus_2;
  PRINTF("det_mod_pk_SE_4x4_invert(): alpha'=%wu %wu | %wu %wu\n",
   invM[10]%p_deg_k,invM[11]%p_deg_k,invM[14]%p_deg_k,invM[15]%p_deg_k);
  PRINTF("det_mod_pk_SE_4x4_invert(): gamma =%wu %wu | %wu %wu\n",
   source_NE_0[0]%p_deg_k,source_NE_0[1]%p_deg_k,source_NE_1[0]%p_deg_k,source_NE_1[1]%p_deg_k);
  det_mod_pk_mul_negate_gamma_alphaINV(
   source_NE_0, source_NE_1,
   invM+10, invM+14,
   mod);
  PRINTF("det_mod_pk_SE_4x4_invert(): -gamma*alpha'=%wu %wu | %wu %wu\n",
   source_NE_0[0]%p_deg_k,source_NE_0[1]%p_deg_k,source_NE_1[0]%p_deg_k,source_NE_1[1]%p_deg_k);
  mp_limb_t* source_SW_0=rows[dim_minus_2  ]+(dim_minus_2-2);
  mp_limb_t* source_SW_1=rows[dim_minus_2+1]+(dim_minus_2-2);
  // SW corner of source: betta -> -betta*epsln'
  det_mod_pk_mul_negate_gamma_alphaINV(
   source_SW_0,source_SW_1,
   invM,invM+4,
   mod);
  // SW corner of target: zetta = alpha' * SW corner of source
  det_mod_pk_mul_2x2(
   invM+8, invM+12,
   invM+10, invM+14,
   source_SW_0,source_SW_1,
   mod);
  // NE corner of target: epsln' * NE corner of source
  det_mod_pk_mul_2x2(
   invM+2, invM+6,
   invM, invM+4,
   source_NE_0, source_NE_1,
   mod);
  // SE corner of target: alpha' -> alpha' + zetta * NE corner of source
  det_mod_pk_mul_add_2x2(
   invM+10, invM+14,
   invM+8, invM+12,
   source_NE_0,source_NE_1,
   mod);
  #undef p_deg_k
 }

mp_limb_t
det_mod_pk_fix_SE_corner(mp_limb_t* I,nmod_mat_t M,slong* negate_det,
 const p_k_pk_t pp)
/*
switch rows so lower-right 4x4 corner is non-singular modulo p,
 or at least put a good element into lower-right corner

return 0 if M is singular
return 1 if the attempt to make the south-east corner good failed
if the attempt succeeded, return 2+determinant(the corner), destroy SE 
 corner of M, set I to inverse of the corner

when switching rows, add 1 to negate_det
*/
 {
  mp_limb_t r,alpha_inv;
  // start with pivot in lower-right corner
  r=I[15]=det_mod_pk_SE_0th_row( M, negate_det, pp.p );
  SHOW_0TH_COL("after 0th row",M);
  if( 0==r )
   {
    I[0]=0; // 0th pivot not invertible
    return det_mod_pk_examine_last_column( M, negate_det, pp );
   }
  // row 1
  if( 0==det_mod_pk_SE_1st_row( I, M, negate_det, pp ) )
   return 1;
  ASSERT(1 == n_mulmod_preinv_4arg(I[0],M->rows[M->r-1][M->r-1], 
   M->mod.n,M->mod.ninv) % pp.p_deg_k);
  SHOW_0TH_COL("after 1st row",M);
  // protect inverse of 0th pivot
  alpha_inv=I[0];
  r=det_mod_pk_SE_row_23( I, M, negate_det, pp.p );
  PRINTF("det_mod_pk_SE_row_23() result: %wu\n",r);
  SHOW_0TH_COL("after row_23",M);
  if( 0==r )
   {
    I[0]=alpha_inv; // 0th pivot inverted
    return 1;
   }
  r=2+det_mod_pk_SE_corner_det( I, M, r );
  PRINTF("det_mod_pk_SE_corner_det(): %wu\n",r-2);
  det_mod_pk_SE_4x4_invert( I, M, pp );
  return r;
 }

static __inline__ void
det_mod_pk_cutoff_1(nmod_mat_t M,mp_limb_t alpha_inv,const p_k_pk_t pp,
  mp_limb_t* scrtch)
/*
divide away degree of p from C

D := D-C*A'*B 
*/
 {
  /*
   for non-static non-inline version of this subroutine length in instructions is
    155+22*2*dim_minus_1 or
    156+22*  dim_minus_1
   (optimized version wins with a score 22*dim_minus_1-1)
  */
  const slong dim_minus_1 = --M->r;
  #if DUMP_cutoff_1_CALL
   flint_printf("%x.",dim_minus_1);
  #endif
  mp_limb_t** rows=M->rows;
  mp_limb_t* B0_ptr=rows[dim_minus_1];
  mp_limb_t t,A;
  mp_limb_t* B_ptr;
  mp_limb_t* D_ptr;
  slong i,j;
  const nmod_t mod=M->mod;
  // Divide by a power of p if necessary
  if( 0 == alpha_inv )
   {
    A=B0_ptr[dim_minus_1];
    t=(mp_limb_t)n_remove( &A, pp.p );
    t=n_pow_speedup(pp.p, t);
    for(j=dim_minus_1;j--;)
     rows[j][dim_minus_1] /= t;
    alpha_inv=inv_mod_pk_3arg(A,pp,mod);
   }
  #if defined(MULADD_pk_pointer)
   //negate scrtch[] and use optimized mul-add instead of unoptimized mul-sub
   A=mod.n-alpha_inv;
   for(j=dim_minus_1;j--;)
    scrtch[j]=n_mulmod_preinv_4arg(rows[j][dim_minus_1],A,mod.n,mod.ninv);
   for(i=dim_minus_1;i--;)
    {
     t=scrtch[i];
     D_ptr=rows[i];
     B_ptr=B0_ptr;
     for(j=dim_minus_1;j--;B_ptr++,D_ptr++)
      MULADD_pk_pointer( D_ptr, t,B_ptr[0], mod.n,mod.ninv);
    }
  #else
   // scrtch[0..dim_minus_1-1] := column M[0..dim_minus_1-1][dim_minus_1]*A
   for(j=dim_minus_1;j--;)
    scrtch[j]=n_mulmod_preinv_4arg(rows[j][dim_minus_1],alpha_inv,
     mod.n,mod.ninv);
   // D := D - scrtch transposed * B
   for(i=dim_minus_1;i--;)
    {
     t=scrtch[i];
     D_ptr=rows[i];
     B_ptr=B0_ptr;
     for(j=dim_minus_1;j--;B_ptr++,D_ptr++)
      D_ptr[0] = n_submod( 
                  D_ptr[0],
                  n_mulmod_preinv_4arg(t,B_ptr[0], mod.n,mod.ninv),
                  mod.n);
    }
  #endif
 }

#if defined(ALIGN_inv_array)
 // WTF? GCC silently ignored the 2 movdqa commands
 #define SWALLOW_sou            \
  __asm__ volatile               \
   (                              \
    "movdqa %2,%0\n\tmovdqa %3,%1" \
    : "=x" (cont0), "=x" (cont1)     \
    : "m" (sou[0]), "m" (sou[2])       \
   );
#else
 #define SWALLOW_sou
#endif

#define A_BY_B( Arow )        \
    tgt=scrtch+Arow;              \
    sou=Ainv+4*Arow;                   \
    SWALLOW_sou                          \
    SCALAR_4(b0[0],b1[0],b2[0],b3[0]);   \
    for(j=1;j<dim_minus_4;j++)           \
     {                                   \
      tgt += 4;                          \
      SCALAR_4(b0[j],b1[j],b2[j],b3[j]); \
     }

/*
SCALAR_4(x0,x1,x2,x3): tgt[0] := x0..x3 DOT sou[0]..sou[3]
 multiple DOT operations for single sou performed
*/

#if defined(ALIGN_inv_array)
#define SCALAR_4(x0,x1,x2,x3 )   \
 {                                            \
  register mp_limb_t t1,t2;                      \
  __asm__                                              \
   (                                                \
    "pextrq $0,%3,%%rax\n\t"                         \
    "mulq %q5\n\t"                                    \
    "movq %%rax,%q2\n\t"                               \
    "movq %%rdx,%q1\n\t"                                \
    "xorq %q0,%q0\n\t"                                   \
    "pextrq $1,%3,%%rax\n\t"                              \
    "mulq %q6\n\t"                                         \
    "addq %%rax,%q2\n\t"                                    \
    "adcq %%rdx,%q1\n\t"                                     \
    "adcq $0x0,%q0\n\t"                                       \
    "pextrq $0,%4,%%rax\n\t"                                   \
    "mulq %q7\n\t"                                              \
    "addq %%rax,%q2\n\t"                                         \
    "adcq %%rdx,%q1\n\t"                                          \
    "adcq $0x0,%q0\n\t"                                           \
    "pextrq $1,%4,%%rax\n\t"                                       \
    "mulq %q8\n\t"                                                  \
    "addq %%rax,%q2\n\t"                                             \
    "adcq %%rdx,%q1\n\t"                                              \
    "adcq $0x0,%q0\n\t"                                                \
    : "=r" (t2), "=r" (t1), "=r" (t)                                   \
    : "x" (cont0), "x" (cont1), "m" (x0), "m" (x1), "m" (x2), "m" (x3) \
    : "rax", "rdx"                                                    \
   );                                                                \
  NMOD_RED3_pk(t2,t1,t, mod.n,mod.ninv,mod.norm);                  \
  tgt[0]=t;                                                     \
 }
#endif

#if VECTOR_DOT_IN_cutoff_4 && !defined(SCALAR_4)
 // This SCALAR_4 and below has same speed
 #define SCALAR_4(b_0,b_1,b_2,b_3)            \
  {                                            \
   VECTOR_DOT_HEAD(sou[0],b_0);                 \
   VECTOR_DOT_BODY(sou[1],b_1);                  \
   VECTOR_DOT_BODY(sou[2],b_2);                   \
   VECTOR_DOT_BODY(sou[3],b_3);                    \
   VECTOR_DOT_TAIL(tgt[0],mod.n,mod.ninv,mod.norm); \
  }
#endif

#if !defined(SCALAR_4)
 #define SCALAR_4(b_0,b_1,b_2,b_3)                     \
  t=n_mulmod_preinv_4arg( sou[0],b_0, mod.n,mod.ninv ); \
  t=n_addmod(                                           \
   t,                                                   \
   n_mulmod_preinv_4arg( sou[1],b_1, mod.n,mod.ninv ),  \
   mod.n);                                              \
  t=n_addmod(                                           \
   t,                                                   \
   n_mulmod_preinv_4arg( sou[2],b_2, mod.n,mod.ninv ),  \
   mod.n);                                              \
  tgt[0]=n_addmod(                                      \
   t,                                                   \
   n_mulmod_preinv_4arg( sou[3],b_3, mod.n,mod.ninv ),  \
   mod.n);
#endif

#define MUL_SUB(r, s, m0,m1)                 \
 r=n_submod(                                  \
  s,                                            \
  n_mulmod_preinv_4arg( m0,m1, mod.n,mod.ninv ), \
  mod.n);

static __inline__ void
det_mod_pk_cutoff_4(mp_limb_t* Ainv,nmod_mat_t M,mp_limb_t* scrtch)
/*
D := D-B*A'*C 
*/
 {
  mp_limb_t** rows=M->rows;
  const slong dim_minus_4=( M->r -= 4);
  const mp_limb_t* b0=rows[dim_minus_4];
  const mp_limb_t* b1=rows[dim_minus_4+1];
  const mp_limb_t* b2=rows[dim_minus_4+2];
  const mp_limb_t* b3=rows[dim_minus_4+3];
  mp_limb_t* tgt;
  mp_limb_t* sou;
  #if defined(ALIGN_inv_array)
   register uint128_t cont0 __asm__("xmm0");
   register uint128_t cont1 __asm__("xmm1");
  #endif
  mp_limb_t* so2;
  const nmod_t mod=M->mod;
  slong i,j;
  mp_limb_t t;
  // scrtch := (A'*B) transposed
  A_BY_B(0);
  A_BY_B(1);
  A_BY_B(2);
  A_BY_B(3);
  // D := D-C * scrtch transposed
  for(i=dim_minus_4;i--;)
   {
    tgt=rows[i];
    sou=tgt+dim_minus_4;
    so2=scrtch;
    /*
     attempt to cache in registers: no gain
      mp_limb_t sou0=sou[0];
      mp_limb_t sou1=sou[1];
      ...
    */
    for(j=dim_minus_4;j--;tgt++,so2 += 4)
     {
      #if defined(VECTOR_DOT_TAIL_add)
       /*
        this should be at shorter, since MUL_SUB is approximately 30
         instructions, VECTOR_DOT_BODY is 5, VECTOR_DOT_TAIL_add is 3+???
         ... hmm, don't know
       */
       VECTOR_DOT_HEAD(sou[0],so2[0]);
       VECTOR_DOT_BODY(sou[1],so2[1]);
       VECTOR_DOT_BODY(sou[2],so2[2]);
       VECTOR_DOT_BODY(sou[3],so2[3]);
       t=mod.n-tgt[0];
       VECTOR_DOT_TAIL_add(t, mod.n,mod.ninv,mod.norm);
       n_negmod_opt(t,mod.n);
       tgt[0]=t;
      #else
       t=tgt[0];
       MUL_SUB(t, t, sou[0],so2[0] );
       MUL_SUB(t, t, sou[1],so2[1] );
       MUL_SUB(t, t, sou[2],so2[2] );
       MUL_SUB(tgt[0], t, sou[3],so2[3] );
      #endif
     }
   }
 }

#undef MUL_SUB
#undef A_BY_B
#undef SCALAR_4

mp_limb_t
nmod_mat_det_mod_pk_4block(nmod_mat_t M,const p_k_pk_t pp,mp_limb_t* scrtch)
 {
  #define p_deg_k pp.p_deg_k
  slong negate_det=0,dim;
  mp_limb_t inv[16]
  #if defined(ALIGN_inv_array)
   __attribute__((aligned(ALIGN_inv_array)))
  #endif
   ;
  mp_limb_t c,result=UWORD_MAX;
  const nmod_t mod=M->mod;
  // Reduce dimension
  while( (dim=M->r) > 4 )
   {
    c=det_mod_pk_fix_SE_corner( inv, M, &negate_det, pp );
    #if BUG_IN_cutoff_4
     flint_printf("dim=%w fix_SE_corner():%wu\n",dim,c);
     show_all_but_SE_4x4_corner(M,p_deg_k);
     show_Ainv(inv,p_deg_k);
    #endif
    if(c==0)
     return 0; // zero column
    if(c==1)
     {
      c=dim-1;
      if(result != UWORD_MAX)
       result=n_mulmod_preinv_4arg( M->rows[c][c], result, mod.n, mod.ninv );
      else
       result=M->rows[c][c];
      assert
       ( 
        ( (0<inv[0]) &&
          (1==n_mulmod_preinv_4arg(M->rows[c][c],inv[0],mod.n,mod.ninv)%p_deg_k)
        )
        || 
        ( (0==inv[0]) &&
          (0==(M->rows[c][c])%pp.p)
        )
       );
      det_mod_pk_cutoff_1( M, inv[0], pp, scrtch );
     }
    else
     {
      if(result != UWORD_MAX)
       result=n_mulmod_preinv_4arg( c-2, result, mod.n, mod.ninv );
      else
       result=c-2;
      det_mod_pk_cutoff_4( inv, M, scrtch );
     }
   }
  // For small dim, use division-less algorithm
  switch(dim)
   {
    case 4 : c = nmod_mat_det_dim4(M); break;
    case 3 : c = nmod_mat_det_dim3(M); break;
    case 2 : c = nmod_mat_det_dim2(M); break;
    default: c = M->rows[0][0];
   }
  #if BUG_IN_cutoff_4
   flint_printf("det_mod_pk(): negate_det=%wu result=%wu c=%wu\n",
    negate_det,result%p_deg_k,c % p_deg_k);
  #endif
  if(negate_det)
   c=n_negmod( c, mod.n );
  if(result != UWORD_MAX)
   {
    c=n_mulmod_preinv_4arg( c, result, mod.n, mod.ninv );
    M->r=M->c;
   }
  if( mod.n & 1 )
   return c;
  return c % p_deg_k;
  #undef p_deg_k
 }

#undef ASSERT
