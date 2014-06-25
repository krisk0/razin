// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// 4x4 determinant calculation is based on comments in Sage file
//  matrix_integer_dense.pyx

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <assert.h>

#define LOUD_4x4_INVERT 0
#if LOUD_4x4_INVERT
 #define PRINTF flint_printf
 #define SHOW_0TH_COL(m,a) flint_printf("%s 0th col: %wu %wu %wu %wu %wu\n", \
  m,nmod_mat_entry(a,0,0),nmod_mat_entry(a,1,0),nmod_mat_entry(a,2,0),      \
  nmod_mat_entry(a,3,0),nmod_mat_entry(a,4,0));
#else
 #define PRINTF
 #define SHOW_0TH_COL
#endif

#define MUL n_mulmod_preinv_4arg
#define MUL_mod_n(x,y) n_mulmod_preinv_4arg(x,y,n,i)
#define SUB_mod_n(x,y) n_submod(x,y,n)
#define ADD_mod_n(x,y) n_addmod(x,y,n)
#define DET_4x4                           \
 a=MUL_mod_n( r0[3],r1[2] );               \
 a=SUB_mod_n( a, MUL_mod_n( r0[2],r1[3] ) );\
 b=MUL_mod_n( r2[1],r3[0] );                 \
 b=SUB_mod_n( b, MUL_mod_n( r2[0],r3[1] ) ); \
 r=MUL_mod_n( a, b );                        \
 a=MUL_mod_n( r0[1],r1[3] );                 \
 a=SUB_mod_n( a, MUL_mod_n( r0[3],r1[1] ) ); \
 b=MUL_mod_n( r2[2],r3[0] );                 \
 b=SUB_mod_n( b, MUL_mod_n( r2[0],r3[2] ) ); \
 r=ADD_mod_n( r, MUL_mod_n(a,b) );           \
 a=MUL_mod_n( r0[2],r1[1] );                 \
 a=SUB_mod_n( a, MUL_mod_n( r0[1],r1[2] ) ); \
 b=MUL_mod_n( r2[3],r3[0] );                 \
 b=SUB_mod_n( b, MUL_mod_n( r2[0],r3[3] ) ); \
 r=ADD_mod_n( r, MUL_mod_n(a,b) );           \
 a=MUL_mod_n( r0[3],r1[0] );                 \
 a=SUB_mod_n( a, MUL_mod_n( r0[0],r1[3] ) ); \
 b=MUL_mod_n( r2[2],r3[1] );                 \
 b=SUB_mod_n( b, MUL_mod_n( r2[1],r3[2] ) ); \
 r=ADD_mod_n( r, MUL_mod_n(a,b) );           \
 a=MUL_mod_n( r0[0],r1[2] );                 \
 a=SUB_mod_n( a, MUL_mod_n( r0[2],r1[0] ) ); \
 b=MUL_mod_n( r2[3],r3[1] );                 \
 b=SUB_mod_n( b, MUL_mod_n( r2[1],r3[3] ) ); \
 r=ADD_mod_n( r, MUL_mod_n(a,b) );           \
 a=MUL_mod_n( r0[1],r1[0] );                 \
 a=SUB_mod_n( a, MUL_mod_n( r0[0],r1[1] ) ); \
 b=MUL_mod_n( r2[3],r3[2] );                 \
 b=SUB_mod_n( b, MUL_mod_n( r2[2],r3[3] ) ); \
 r=ADD_mod_n( r, MUL_mod_n(a,b) );

mp_limb_t
nmod_mat_det_dim4(nmod_mat_t M)
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
  mp_limb_t** rows=M->rows;
  mp_limb_t* r0=rows[0];
  mp_limb_t* r1=rows[1];
  mp_limb_t* r2=rows[2];
  mp_limb_t* r3=rows[3];
  mp_limb_t n=M->mod.n;
  mp_limb_t i=M->mod.ninv;
  mp_limb_t r,a,b;
  DET_4x4
  return r;
 }

mp_limb_t
nmod_mat_det_dim4_SE(nmod_mat_t M)
/*
count determinant of 4x4 submatrice in lower-right corner
use same algorithm as nmod_mat_det_dim4()
*/
 {
  slong s=M->r-4;
  mp_limb_t** rows=M->rows;
  mp_limb_t* r0=rows[s]+s;
  mp_limb_t* r1=rows[s+1]+s;
  mp_limb_t* r2=rows[s+2]+s;
  mp_limb_t* r3=rows[s+3]+s;
  mp_limb_t n=M->mod.n;
  mp_limb_t i=M->mod.ninv;
  mp_limb_t r,a,b;
  DET_4x4
  return r;
 }

#define DIM2_DET( rez, n, i, r0, r1 )     \
  rez=          MUL( r0[0], r1[1], n, i ); \
  rez=n_submod( rez,                       \
                MUL( r0[1], r1[0], n, i ), \
            n );

mp_limb_t
nmod_mat_det_dim3(nmod_mat_t A)
 {
  mp_limb_t** rows=A->rows;
  mp_limb_t* r0=rows[0];
  mp_limb_t* r1=rows[1];
  mp_limb_t* r2=rows[2];
  mp_limb_t n=A->mod.n;
  mp_limb_t i=A->mod.ninv;
  mp_limb_t rez,t;
  DIM2_DET( rez, n, i, r0, r1 );
  rez=MUL( rez, r2[2], n, i );
  DIM2_DET( t, n, i, r0, r2 );
  t=MUL( t, r1[2], n, i );
  rez=n_submod( rez, t, n);
  DIM2_DET( t, n, i, r1, r2 );
  t=MUL( t, r0[2], n, i );
  rez=n_addmod( rez, t, n);
  return rez;
 }

mp_limb_t
nmod_mat_det_dim2(nmod_mat_t A)
 {
  mp_limb_t* r0=A->rows[0];
  mp_limb_t* r1=A->rows[1];
  mp_limb_t n=A->mod.n;
  mp_limb_t i=A->mod.ninv;
  mp_limb_t rez;
  DIM2_DET( rez, n, i, r0, r1 );
  return rez;
 }

#undef DIM2_DET
#undef MUL
#undef MUL_mod_n
#undef ADD_mod_n
#undef SUB_mod_n
#undef DET_4x4

mp_limb_t det_mod_pk_cutoff_1(nmod_mat_t M,mp_limb_t p,mp_limb_t p_deg_k,
 ulong k,mp_limb_t* scrtch);
mp_limb_t det_mod_pk_cutoff_4(mp_limb_t* I,nmod_mat_t M,mp_limb_t SE_det,
 mp_limb_t p,mp_limb_t p_deg_k,ulong k,mp_limb_t* scrtch);

static __inline__ mp_limb_t
det_mod_pk_SE_0th_row(nmod_mat_t M,slong* negate_det,mp_limb_t p,
  mp_limb_t p_deg_k)
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
  tempp=n_mulmod_preinv_4arg(tempp, zetta, mod.n, mod.ninv);
  invM[15]=n_addmod( tempp, alpha_inv, mod.n );
  #if LOUD_4x4_INVERT
   flint_printf("SE 2x2 corner inverse: %wu %wu | %wu %wu\n",
    invM[10],invM[11],invM[14],invM[15]);
  #endif
 }

static __inline__ mp_limb_t 
det_mod_pk_SE_1st_row(mp_limb_t* invM,nmod_mat_t M,slong* negate_det,
  mp_limb_t p,ulong k,mp_limb_t p_deg_k)
// return 1 on success  
 {
  mp_limb_t** rows=M->rows;
  const mp_limb_t dim_minus_1=M->r-1;
  const mp_limb_t dim_minus_2=dim_minus_1-1;
  const nmod_t mod=M->mod;
  const mp_limb_t alpha=invM[15];
  const mp_limb_t alpha_inv=inv_mod_pk(alpha,p,k,p_deg_k,mod.n,mod.ninv);
  const mp_limb_t betta=rows[dim_minus_1][dim_minus_2];
  const mp_limb_t alpha_inv_by_beta=n_mulmod_preinv_4arg(alpha_inv,betta,
   mod.n, mod.ninv );
  mp_limb_t* tail;
  mp_limb_t gamma,delta,epsln;
  assert( n_mulmod_preinv_4arg(alpha,alpha_inv,mod.n,mod.ninv) % p_deg_k==1 );
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
    if(epsln % p)
     {
      if(i != dim_minus_2)
       {
        (*negate_det) ^= 1;
        MP_PTR_SWAP( rows[i], rows[dim_minus_2] );
       }
      invM[12]=epsln; // put the pivot into SW corner of invM
      epsln=inv_mod_pk(epsln,p,k,p_deg_k,mod.n,mod.ninv);
      det_mod_pk_SE_2x2_invert( invM, M, alpha_inv, betta, epsln, gamma );
      return 1;
     }
   }
  return 0;
 }

#define MUL(x,y) n_mulmod_preinv_4arg(x,y,mod.n,mod.ninv)
#define ADD(x,y) n_addmod(x,y,mod.n)
#define SUB(x,y) n_submod(x,y,mod.n)

static void 
det_mod_pk_mul_2x2( 
  mp_limb_t* r0,mp_limb_t* r1,
  const mp_limb_t* a0,const mp_limb_t* a1,
  const mp_limb_t* b0,const mp_limb_t* b1,
  const nmod_t mod)
 {
  r0[0]=ADD( MUL(a0[0],b0[0]), MUL(a0[1],b1[0]) );
  r0[1]=ADD( MUL(a0[0],b0[1]), MUL(a0[1],b1[1]) );
  r1[0]=ADD( MUL(a1[0],b0[0]), MUL(a1[1],b1[0]) );
  r1[1]=ADD( MUL(a1[0],b0[1]), MUL(a1[1],b1[1]) );
 }

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
  eps0[0]=SUB( sou0[0], ADD( MUL(sou0[2],zet0[0]), MUL(sou0[3],zet1[0]) ) );
  eps0[1]=SUB( sou0[1], ADD( MUL(sou0[2],zet0[1]), MUL(sou0[3],zet1[1]) ) );
  eps1[0]=SUB( sou1[0], ADD( MUL(sou1[2],zet0[0]), MUL(sou1[3],zet1[0]) ) );
  eps1[1]=SUB( sou1[1], ADD( MUL(sou1[2],zet0[1]), MUL(sou1[3],zet1[1]) ) );
 }

#define ROW_23( i, j )                                          \
 cand_0=rows[i]+dim_minus_4;                                          \
 cand_1=rows[j]+dim_minus_4;                                              \
 det_mod_pk_mul_sub_2x2(invM,epsi1,   cand_0,cand_1,   zeta0,zeta1,   mod); \
 PRINTF("i=%w j=%w epsi=%wu %wu | %wu %wu\n",i,j,       \
  invM[0],invM[1],epsi1[0],epsi1[1]);                   \
 ok=SUB( MUL(invM[0],epsi1[1]), MUL(invM[1],epsi1[0]) );                     \
 PRINTF("epsi det modulo %wu = %wu\n",p,ok%p);           \
 if(0 == ok % p )                                                             \
  ok=0;

static __inline__ mp_limb_t 
det_mod_pk_SE_row_23(mp_limb_t* invM,nmod_mat_t M,slong* negate_det,
  mp_limb_t p,ulong k,mp_limb_t p_deg_k)
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

static __inline__ mp_limb_t
det_mod_pk_examine_last_column(nmod_mat_t M,slong* negate_det,mp_limb_t p,
  mp_limb_t p_deg_k)
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
    e=rows[i][shi] % p_deg_k;
    if(e)
     {
      j=n_remove( &e, p );// TODO: maybe use n_remove2_precomp() here?
      if(j<all_zero)
       {
        all_zero=j;
        best_row=i;
        if( 1==j )
         break;
       }
     }
   }
  if(best_row != shi)
   {
    (*negate_det) ^= 1;
    MP_PTR_SWAP( rows[best_row], rows[shi] );
   }
  return (best_row != -1);
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
nmod_invert_2x2_5arg(mp_limb_t* s0,mp_limb_t* s1,const nmod_t mod,
  mp_limb_t p,ulong k,mp_limb_t p_deg_k)
 {
  mp_limb_t delta,gamma,betta,alpha=s1[1];
  mp_limb_t epsil,zetta,theta; 
  if(alpha % p)
   {
    betta=s1[0];
    delta=s0[0];
    gamma=s0[1];
    alpha=inv_mod_pk(alpha,p,k,p_deg_k,mod.n,mod.ninv);
    // epsil=delta-gamma*alpha'*betta
    epsil=n_mulmod_preinv_4arg(gamma,alpha,mod.n,mod.ninv);
    epsil=n_mulmod_preinv_4arg(epsil,betta,mod.n,mod.ninv);
    epsil=n_submod(delta,epsil,mod.n);
    // zetta=-epsil'*alpha'*betta
    s0[0]=epsil=inv_mod_pk(epsil,p,k,p_deg_k,mod.n,mod.ninv);
    assert(epsil<p_deg_k);
    zetta=n_mulmod_preinv_4arg(p_deg_k-epsil,alpha,mod.n,mod.ninv);
    s1[0]=zetta=n_mulmod_preinv_4arg(zetta,betta,mod.n,mod.ninv);
    // theta=gamma*alpha'
    theta=n_mulmod_preinv_4arg(p_deg_k-alpha,gamma,mod.n,mod.ninv);
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
    gamma=inv_mod_pk(gamma,p,k,p_deg_k,mod.n,mod.ninv);
    // epsil=delta-gamma*alpha'*betta
    epsil=n_mulmod_preinv_4arg(gamma,alpha,mod.n,mod.ninv);
    epsil=n_mulmod_preinv_4arg(epsil,betta,mod.n,mod.ninv);
    epsil=n_submod(delta,epsil,mod.n);
    s0[1]=epsil=inv_mod_pk(epsil,p,k,p_deg_k,mod.n,mod.ninv);
    assert(epsil<p_deg_k);
    zetta=n_mulmod_preinv_4arg(p_deg_k-epsil,gamma,mod.n,mod.ninv);
    s1[1]=zetta=n_mulmod_preinv_4arg(zetta,betta,mod.n,mod.ninv);
    theta=n_mulmod_preinv_4arg(p_deg_k-gamma,alpha,mod.n,mod.ninv);
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
 }

void
det_mod_pk_mul_add_2x2(
  mp_limb_t* tgt0,mp_limb_t* tgt1,
  mp_limb_t* rho0,mp_limb_t* rho1,
  mp_limb_t* eta0,mp_limb_t* eta1,
  const nmod_t mod)
// tgt := tgt + rho * eta
 {
  tgt0[0] = n_addmod( tgt0[0],
    n_mulmod_preinv_4arg(rho0[0],eta0[0], mod.n, mod.ninv),
    mod.n
   );
  tgt0[0] = n_addmod( tgt0[0],
    n_mulmod_preinv_4arg(rho0[1],eta1[0], mod.n, mod.ninv),
    mod.n
   );
  tgt0[1] = n_addmod( tgt0[1],
    n_mulmod_preinv_4arg(rho0[0],eta0[1], mod.n, mod.ninv),
    mod.n
   );
  tgt0[1] = n_addmod( tgt0[1],
    n_mulmod_preinv_4arg(rho0[1],eta1[1], mod.n, mod.ninv),
    mod.n
   );
  tgt1[0] = n_addmod( tgt1[0],
    n_mulmod_preinv_4arg(rho1[0],eta0[0], mod.n, mod.ninv),
    mod.n
   );
  tgt1[0] = n_addmod( tgt1[0],
    n_mulmod_preinv_4arg(rho1[1],eta1[0], mod.n, mod.ninv),
    mod.n
   );
  tgt1[1] = n_addmod( tgt1[1],
    n_mulmod_preinv_4arg(rho1[0],eta0[1], mod.n, mod.ninv),
    mod.n
   );
  tgt1[1] = n_addmod( tgt1[1],
    n_mulmod_preinv_4arg(rho1[1],eta1[1], mod.n, mod.ninv),
    mod.n
   );
 }

#define SHOW_VECTOR_4(i)                               \
 v=a->rows[i]+(a->r-4);                                 \
 flint_printf("%wu %wu %wu %wu\n",v[0]%n,v[1]%n,v[2]%n,v[3]%n);

void show_SE_4x4_corner(nmod_mat_t a,mp_limb_t n)
 {
  mp_limb_t* v;
  SHOW_VECTOR_4(a->r-4)
  SHOW_VECTOR_4(a->r-3)
  SHOW_VECTOR_4(a->r-2)
  SHOW_VECTOR_4(a->r-1)
 }

static __inline__ void
det_mod_pk_SE_4x4_invert(mp_limb_t* invM,nmod_mat_t M,mp_limb_t p,ulong k,
  mp_limb_t p_deg_k)
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
  nmod_invert_2x2_5arg(invM,invM+4,mod,p,k,p_deg_k);
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
 }

mp_limb_t
det_mod_pk_fix_SE_corner(mp_limb_t* I,nmod_mat_t M,slong* negate_det,
 mp_limb_t p,ulong k,mp_limb_t p_deg_k)
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
  mp_limb_t r;
  // start with pivot in lower-right corner
  r=I[15]=det_mod_pk_SE_0th_row( M, negate_det, p, p_deg_k );
  SHOW_0TH_COL("after 0th row",M);
  if( 0==r )
   return det_mod_pk_examine_last_column( M, negate_det, p, p_deg_k );
  assert( r == M->rows[M->r-1][M->r-1] );
  // row 1
  if( 0==det_mod_pk_SE_1st_row( I, M, negate_det, p, k, p_deg_k ) )
   return 1;
  SHOW_0TH_COL("after 1st row",M);
  assert( M->rows[M->r-1][M->r-1] );
  r=det_mod_pk_SE_row_23( I, M, negate_det, p, k, p_deg_k );
  PRINTF("det_mod_pk_SE_row_23() result: %wu\n",r);
  SHOW_0TH_COL("after row_23",M);
  if( 0==r )
   return 1;
  assert( M->rows[M->r-1][M->r-1] );
  r=2+det_mod_pk_SE_corner_det( I, M, r );
  PRINTF("det_mod_pk_SE_corner_det(): %wu\n",r-2);
  det_mod_pk_SE_4x4_invert( I, M, p, k, p_deg_k );
  return r;
 }

mp_limb_t
nmod_mat_det_mod_pk(nmod_mat_t M,mp_limb_t p,mp_limb_t p_deg_k,ulong k,
                    mp_limb_t* scrtch)
 {
  slong negate_det=0;
  slong dim=M->r;
  mp_limb_t rez;
  if(dim <= 4)
   {
    // Use division-less algorithm
    if(dim == 4)
     rez = nmod_mat_det_dim4(M);
    else
     {
      if(dim == 3)
       rez = nmod_mat_det_dim3(M);
      else
       {
        if(dim == 2)
         rez = nmod_mat_det_dim2(M);
        else
         rez = M->entries[0];
       }
     }
    return rez % p_deg_k;
   }
  mp_limb_t inv[16];
  rez=det_mod_pk_fix_SE_corner( inv, M, &negate_det, p, k, p_deg_k );
  if(rez==0)
   return rez;
  printf("end of rails\n");
  abort();
  if(rez==1)
   rez=det_mod_pk_cutoff_1(M,p,p_deg_k,k,scrtch);
  else
   rez=det_mod_pk_cutoff_4(inv,M,rez-2,p,p_deg_k,k,scrtch);
  if(negate_det)
   rez=n_negmod(rez,p_deg_k);
  return rez;
 }
