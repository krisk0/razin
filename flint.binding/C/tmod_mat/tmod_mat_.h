// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// t_invmod() macro uses a piece of public domain code owned 
//  by Jean-Guillaume Dumas

#ifndef TMOD_MAT__H
#define TMOD_MAT__H

#include <flint/flint.h>
#include <flint/fmpz_mat.h>

typedef struct
 {
  slong r;
  slong c;
  mp_limb_t**   rows;
  mp_limb_t*  entries;
 }
tmod_mat_struct;
typedef tmod_mat_struct tmod_mat_t[1];
#define tmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])

void tmod_mat_init_fast(tmod_mat_t mat, long rows, long cols);
void tmod_mat_clear(tmod_mat_t mat);
void tmod_mat_print_hex(char* m,const tmod_mat_t S);
mp_limb_t dot_modulo_t_kind0(mp_limb_t* a, fmpz* b, slong n);
void tmod_mat_init_set_fmpz_mat(tmod_mat_t tgt, const fmpz_mat_t sou);
mp_limb_t det_mod_t(const fmpz_mat_t S);

// TODO: this number-oriented subroutine should be in /ulong_extras
static mp_limb_t inv_mod__2_64__tab[13]={
 1,0xAAAAAAAAAAAAAAAB,0xCCCCCCCCCCCCCCCD,0x6DB6DB6DB6DB6DB7,0x8E38E38E38E38E39
 ,0x2E8BA2E8BA2E8BA3,0x4EC4EC4EC4EC4EC5,0xEEEEEEEEEEEEEEEF,0xF0F0F0F0F0F0F0F1,
 0x86BCA1AF286BCA1B,0xCF3CF3CF3CF3CF3D,0xD37A6F4DE9BD37A7,0x8F5C28F5C28F5C29};

static __inline__ mp_limb_t 
t_invmod(mp_limb_t a)
 {
  #define m 64
  mp_limb_t s;
  mp_limb_t t=a-1;
  if(t<26)
   return inv_mod__2_64__tab[ t>>1 ];
  // inspired by public domain code created by Jean-Guillaume Dumas
  count_trailing_zeros(s,t);
  t >>= s; 
  //assert( a == (1<<s)*t + 1 );
  //assert( t&1 );
  mp_limb_t U=2-a;
  mp_limb_t amone=a-1;
  long i,i_max;
  i_max=m/s;
  for(i=1;i<i_max;i <<= 1)
   {
    amone *= amone;
    U *= amone+1;
   }
  #undef m
  return U;
 }

static __inline__ void
tmod_vec_mul(mp_limb_t* tgt,mp_limb_t* sou,slong siz,mp_limb_t quo)
 {
  for(;siz--;)
   tgt[siz]=sou[siz]*quo;
 }

static __inline__ void
tmod_vec_mul_3arg(mp_limb_t* vec,slong siz,mp_limb_t quo)
 {
  for(;siz--;)
   vec[siz] *= quo;
 }

static __inline__ void
tmod_vec_neg_2arg(mp_limb_t* vec,slong siz)
 {
  for(;siz--;)
   vec[siz] = -vec[siz];
 }

static __inline__ void
tmod_vec_muladd(mp_limb_t* tgt,mp_limb_t* sou,slong siz,mp_limb_t quo)
 {
  for(;siz--;)
   tgt[siz] += sou[siz]*quo;
 }

#endif
