// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include "../ulong_extras/longlong_.h"

void nmod_mat_mod_t_half(nmod_mat_t tgt, fmpz_mat_t sou);

#if (GMP_LIMB_BITS == 64 && defined (__amd64__))

 #define VECTOR_DOT_HEAD(alpha, betta)  \
   register mp_limb_t Vhi=0,Vmi,Vlo,V1,V0;  \
   umul_ppmm( Vmi,Vlo, alpha,betta );        \

 #define VECTOR_DOT_INIT                 \
   register mp_limb_t Vhi,Vmi,Vlo,V1,V0; 

 #define VECTOR_DOT_STRT(alpha, betta)  \
   Vhi=0;                              \
   umul_ppmm( Vmi,Vlo, alpha,betta );

 #define VECTOR_DOT_BODY(alpha, betta)  \
   umul_ppmm( V1,V0, alpha,betta );     \
   add_sssaa( Vhi,Vmi,Vlo, V1,V0 );     
 
 #define VECTOR_DOT_TAIL(rez, n,ninv,two_128_mod_n)  \
   NMOD_RED3_pk( Vhi,Vmi,Vlo, n,ninv,two_128_mod_n );  \
   rez=Vlo;                                             

 #define VECTOR_DOT_TAIL_easy(rez, n,ninv)                \
   NMOD_RED3_pk_easy( Vhi,Vmi,Vlo, n,ninv );               \
   rez=Vlo;                                                

 // r=s-dot result
 #define VECTOR_DOT_TAIL_sub(r, s, n,ninv,two_128_mod_n)  \
  {                                                      \
   NMOD_RED3_pk( Vhi,Vmi,Vlo, n,ninv,two_128_mod_n );   \
   r=n_submod( s,Vlo, n );                             \
  }
 
 #define VECTOR_DOT_TAIL_sub_easy(r, s, n,ninv)           \
  {                                                      \
   NMOD_RED3_pk_easy( Vhi,Vmi,Vlo, n,ninv          );   \
   r=n_submod( s,Vlo, n );                             \
  }

 // r=r+dot result, if r+dot result < 2**129
 #define VECTOR_DOT_TAIL_add(r, n,ninv)     \
  {                                           \
   add_sssa(Vhi,Vmi,Vlo, r);                    \
   NMOD_RED3_pk_easy( Vhi,Vmi,Vlo, n,ninv );     \
   r=Vlo;                                       \
  }
 
 /*
  VECTOR_DOT_TAIL_3arg() is slower than VECTOR_DOT_TAIL
  4x4 by 4x100 matrix multiplication slowdown: 17%
 */
 #define VECTOR_DOT_TAIL_3arg(rez, n,ninv)             \
  {                                                   \
   NMOD_RED3_pk_5arg( Vhi,Vmi,Vlo, n,ninv );         \
   rez=Vlo;                                         \
  }

#endif
