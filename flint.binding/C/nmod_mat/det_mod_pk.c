// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// 4x4 determinant calculation is based on comments in Sage file
//  matrix_integer_dense.pyx

#include <flint/flint.h>
#include <flint/nmod_mat.h>

#define MUL n_mulmod_preinv_4arg
#define MUL_mod_n(x,y) n_mulmod_preinv_4arg(x,y,n,i)
#define SUB_mod_n(x,y) n_submod(x,y,n)
#define ADD_mod_n(x,y) n_addmod(x,y,n)

mp_limb_t nmod_mat_det_dim4(nmod_mat_t M)
/*
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
  //a = x[3]*x[6]; a-=x[2]*x[7]; b = x[ 9]*x[12]; b-=x[8]*x[13];  r  = a*b
  a=MUL_mod_n( r0[3],r1[2] );
  a=SUB_mod_n( a, MUL_mod_n( r0[2],r1[3] ) );
  b=MUL_mod_n( r2[1],r3[0] );
  b=SUB_mod_n( b, MUL_mod_n( r2[0],r3[1] ) );
  r=MUL_mod_n( a, b );
  //a = x[1]*x[7]; a-=x[3]*x[5]; b = x[10]*x[12]; b-=x[8]*x[14];  r += a*b
  a=MUL_mod_n( r0[1],r1[3] );
  a=SUB_mod_n( a, MUL_mod_n( r0[3],r1[1] ) );
  b=MUL_mod_n( r2[2],r3[0] );
  b=SUB_mod_n( b, MUL_mod_n( r2[0],r3[2] ) );
  r=ADD_mod_n( r, MUL_mod_n(a,b) );
  //a = x[2]*x[5]; a-=x[1]*x[6]; b = x[11]*x[12]; b-=x[8]*x[15];  r += a*b
  a=MUL_mod_n( r0[2],r1[1] );
  a=SUB_mod_n( a, MUL_mod_n( r0[1],r1[2] ) );
  b=MUL_mod_n( r2[3],r3[0] );
  b=SUB_mod_n( b, MUL_mod_n( r2[0],r3[3] ) );
  r=ADD_mod_n( r, MUL_mod_n(a,b) );
  //a = x[3]*x[4]; a-=x[0]*x[7]; b = x[10]*x[13]; b-=x[9]*x[14];  r += a*b
  a=MUL_mod_n( r0[3],r1[0] );
  a=SUB_mod_n( a, MUL_mod_n( r0[0],r1[3] ) );
  b=MUL_mod_n( r2[2],r3[1] );
  b=SUB_mod_n( b, MUL_mod_n( r2[1],r3[2] ) );
  r=ADD_mod_n( r, MUL_mod_n(a,b) );
  //a = x[0]*x[6]; a-=x[2]*x[4]; b = x[11]*x[13]; b-=x[9]*x[15];  r += a*b
  a=MUL_mod_n( r0[0],r1[2] );
  a=SUB_mod_n( a, MUL_mod_n( r0[2],r1[0] ) );
  b=MUL_mod_n( r2[3],r3[1] );
  b=SUB_mod_n( b, MUL_mod_n( r2[1],r3[3] ) );
  r=ADD_mod_n( r, MUL_mod_n(a,b) );
  //a = x[1]*x[4]; a-=x[0]*x[5]; b = x[11]*x[14]; b-=x[10]*x[15]; r += a*b
  a=MUL_mod_n( r0[1],r1[0] );
  a=SUB_mod_n( a, MUL_mod_n( r0[0],r1[1] ) );
  b=MUL_mod_n( r2[3],r3[2] );
  b=SUB_mod_n( b, MUL_mod_n( r2[2],r3[3] ) );
  r=ADD_mod_n( r, MUL_mod_n(a,b) );
  return r;
 }

#define DIM2_DET( rez, n, i, r0, r1 )     \
  rez=          MUL( r0[0], r1[1], n, i ); \
  rez=n_submod( rez,                       \
                MUL( r0[1], r1[0], n, i ), \
            n );

mp_limb_t nmod_mat_det_dim3(nmod_mat_t A)
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

mp_limb_t nmod_mat_det_dim2(nmod_mat_t A)
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
    if(dim==4)
     rez = nmod_mat_det_dim4(M);
    else
     {
      if(dim==3)
       rez = nmod_mat_det_dim3(M);
      else
       {
        if(dim==2)
         rez = nmod_mat_det_dim2(M);
        else
         rez = M->entries[0];
       }
     }
    return rez % p_deg_k;
   }
  rez=0;
  if(negate_det)
   rez=n_negmod(rez,p_deg_k);
  return rez;
 }
