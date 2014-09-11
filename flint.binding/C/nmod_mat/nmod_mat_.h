// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef NMOD_MAT__H
#define NMOD_MAT__H

#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include "../ulong_extras/ulong_extras_.h"

void nmod_mat_mod_t_half(nmod_mat_t tgt, const fmpz_mat_t sou);
void nmod_mat_init_3arg(nmod_mat_t mat, slong r, slong c);
void nmod_mat_transpose_square_tgt_virgin(nmod_mat_t tgt,const nmod_mat_t sou);
void nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim);
mp_limb_t nmod_mat_det_mod_pk_4block(nmod_mat_t M,const p_k_pk_t pp,
 mp_limb_t* scrtch);
void nmod_mat_mul_vec_left(mp_limb_t* r, const mp_limb_t* s, const nmod_mat_t m);
void nmod_mat_mul_pk_classical(nmod_mat_t C,nmod_mat_t A,nmod_mat_t B);
int nmod_mat_is_one(const nmod_mat_t m);
void nmod_mat_copy_entries(nmod_mat_t r,const nmod_mat_t s);
slong nmod_mat_HNF(nmod_mat_t r);
mp_limb_t nmod_mat_diag_product_ZZ_ui(const nmod_mat_t m);

#if GMP_LIMB_BITS == 64 && defined (__amd64__)

 // allowing V2_me/V2_lo to be memory results in invalid code and crash
 //  under gcc 4.7.4 and 4.8.3
 #define VECTOR_DOT_2(tgt, x0,y0, x1,y1, mod)   \
  {                                                    \
   register mp_limb_t V2_dx asm ("rdx");                    \
   register mp_limb_t V2_ax asm ("rax");                       \
   mp_limb_t V2_me,V2_lo; /* 1 register for n + theese 2 = 3 */  \
   register mp_limb_t n=mod.n;                                    \
   asm                                                            \
    (                                                             \
        "mov  %q4,%q1\n\t"                                        \
        "mulq %q5\n\t"                                            \
        "mov  %q1,%q3\n\t"                                        \
        "mov  %q0,%q2\n\t"                                        \
        "mov  %q6,%q1\n\t"                                        \
        "mulq %q7\n\t"                                            \
        "add  %q1,%q3\n\t"                                        \
        "adc  %q0,%q2\n\t"                                        \
        "jnc  1f\n\t"                                             \
        /* maybe subtract n from V2_me */                         \
        "xor  %q1,%q1\n\t"                                        \
        "cmp  %q2,%q8\n\t"                                        \
        "cmovb %q8,%q1\n\t"   /* if V2_me >= n then ax=n */       \
        "sub  %q1,%q2\n\t"                                        \
        /* always subtract n from V2_me */                        \
        "sub  %q8,%q2\n\t"                                        \
     "1: xor  %q0,%q0\n\t"                                        \
        "cmp  %q2,%q8\n\t"                                        \
        "cmovb %q8,%q0\n\t"   /* if V2_me >= n then dx=n */       \
        "sub %q0,%q2\n\t"                                         \
     : "=&d" (V2_dx), "=&a" (V2_ax), "=&r" (V2_me), "=&r" (V2_lo) \
     : "m" (x0),"m" (y0), "m" (x1),"m" (y1), "r" (n)             \
    );                                                         \
   NMOD_RED2_pk_4arg(V2_me,V2_lo, n,mod.ninv);               \
   tgt=V2_lo;                                             \
  }

 // tgt += x0*y0+x1*y1, tgt is memory, 2 instructions longer than 
 //  VECTOR_DOT_2()
 #define VECTOR_DOT_2_add(tgt, x0,y0, x1,y1, mod)  \
  {                                                     \
   register mp_limb_t V2_dx asm ("rdx");                    \
   register mp_limb_t V2_ax asm ("rax");                       \
   mp_limb_t V2_me,V2_lo; /* 1 register for n + theese 2 = 3 */  \
   register mp_limb_t n=mod.n;                                    \
   asm                                                            \
    (                                                             \
        "mov  %q4,%q1\n\t"                                        \
        "mulq %q5\n\t"                                            \
        "add  %q9,%q1\n\t"                                        \
        "adc  $0,%q0\n\t"                                         \
        "mov  %q1,%q3\n\t"                                        \
        "mov  %q0,%q2\n\t"                                        \
        "mov  %q6,%q1\n\t"                                        \
        "mulq %q7\n\t"                                            \
        "add  %q1,%q3\n\t"                                        \
        "adc  %q0,%q2\n\t"                                        \
        "jnc  1f\n\t"                                             \
        /* maybe subtract n from V2_me */                         \
        "xor  %q1,%q1\n\t"                                        \
        "cmp  %q2,%q8\n\t"                                        \
        "cmovb %q8,%q1\n\t"   /* if V2_me >= n then ax=n */       \
        "sub  %q1,%q2\n\t"                                        \
        /* always subtract n from V2_me */                        \
        "sub  %q8,%q2\n\t"                                        \
     "1: xor  %q0,%q0\n\t"                                        \
        "cmp  %q2,%q8\n\t"                                        \
        "cmovb %q8,%q0\n\t"   /* if V2_me >= n then dx=n */       \
        "sub %q0,%q2\n\t"                                         \
     : "=&d" (V2_dx), "=&a" (V2_ax), "=&r" (V2_me), "=&r" (V2_lo) \
     : "m" (x0),"m" (y0), "m" (x1),"m" (y1), "r" (n), "m" (tgt)  \
    );                                                         \
   NMOD_RED2_pk_4arg(V2_me,V2_lo, n,mod.ninv);               \
   tgt=V2_lo;                                             \
  }

 #define VECTOR_DOT_HEAD(alpha, betta)  \
   register mp_limb_t Vhi=0,Vmi,Vlo;     \
   umul_ppmm( Vmi,Vlo, alpha,betta );        

 #define VECTOR_DOT_BODY(alpha, betta) \
   asm                                  \
    (                                   \
     "movq %q3,%%rax\n\t"                \
     "mulq %q4\n\t"                      \
     "addq %%rax,%q2\n\t"                  \
     "adcq %%rdx,%q1\n\t"                   \
     "adcq $0x0,%q0"                        \
     : "+rm" (Vhi), "+rm" (Vmi), "+rm" (Vlo) \
     : "rm" (alpha), "rm" (betta)             \
     : "rax", "rdx"                          \
    ); 
  
 #define VECTOR_DOT_HEAD_greedy(alpha, betta)  \
   register mp_limb_t Vhi=0,Vmi,Vlo,V1,V0;    \
   umul_ppmm( Vmi,Vlo, alpha,betta );        \
 
 #define VECTOR_DOT_BODY_greedy(alpha, betta)  \
   {                                     \
    umul_ppmm( V1,V0, alpha,betta );     \
    add_sssaa( Vhi,Vmi,Vlo, V1,V0 );     \
   }
  
 #define VECTOR_DOT_TAIL(rez, n,ninv,two_128_mod_n)  \
  {                                                   \
   NMOD_RED3_pk( Vhi,Vmi,Vlo, n,ninv,two_128_mod_n );  \
   rez=Vlo;                                             \
  }

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

 // r=r+dot result
 #define VECTOR_DOT_TAIL_add(r, n,ninv,two_128_mod_n) \
  {                                                     \
   add_sssa(Vhi,Vmi,Vlo, r);                             \
   NMOD_RED3_pk( Vhi,Vmi,Vlo, n,ninv,two_128_mod_n );     \
   r=Vlo;                                                \
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

 #define VECTOR_DOT_HEAD_tiny(alpha,betta, n,ninv) \
  mp_limb_t Vt=n_mulmod_preinv_4arg(alpha,betta,n,ninv);
 
 #define VECTOR_DOT_BODY_tiny(alpha,betta, n,ninv) \
  MULADD_pk(Vt, alpha,betta, n,ninv);

 #define VECTOR_DOT_TAIL_tiny(rez) \
  rez=Vt;
  
#endif

#endif
