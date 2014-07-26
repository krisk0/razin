// This program is part of RAZIN
/*
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009, 2010, 2012 William Hart
   Copyright 2011 Fredrik Johansson
*/
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL) version 2.1 or later

// This file contains a modified fragment of FLINT longlong.h, owned by W.Hart 
//  and F.Johansson, and a modified fragment of FLINT mulmod_preinv.c, owned by
//  W.Hart

#ifndef LONGLONG__H
#define LONGLONG__H

#include <flint/flint.h>
#include <flint/nmod_vec.h>

#if (GMP_LIMB_BITS == 64 && defined (__amd64__))

// slim version of add_sssaaaaaa macro found in FLINT longlong.h
#define add_sssaaa0aa(sh, sm, sl, ah, am, al,    bm, bl)  \
  __asm__ ("addq %8,%q2\n\tadcq %6,%q1\n\tadcq %4,%q0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "N"   (0x0),              \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  

#define     add_sssaa(sh, sm, sl,                bm, bl)  \
  __asm__ ("addq %4,%q2\n\tadcq %3,%q1\n\tadcq $0x0,%q0"    \
       : "+r" (sh), "+r" (sm), "+r" (sl)                      \
       : "rme" ((mp_limb_t)(bm)), "rme" ((mp_limb_t)(bl))) 

#define     add_sssa(sh, sm, sl,                    bl)  \
  __asm__ ("addq %3,%q2\n\tadcq $0x0,%q1\n\tadcq $0x0,%q0"    \
       : "+r" (sh), "+r" (sm), "+r" (sl)                      \
       :                          "rme" ((mp_limb_t)(bl))) 

static __inline__
mp_limb_t NMOD_RED2_pk_func(mp_limb_t p,mp_limb_t r,mp_limb_t n,mp_limb_t ninv)
// this function is a piece of n_mulmod_preinv_4arg() which is a modification
//  of FLINT n_mulmod_preinv()
 {
  mp_limb_t q1,q0;
  umul_ppmm(q1, q0, ninv, p);
  add_ssaaaa(q1, q0, q1, q0, p, r);
  r -= (q1 + 1) * n;
  if (r >= q0)
   r += n;
  return (r < n ? r : r - n);
 }

// result in r, q*: scratch
#define NMOD_RED2_pk( p,r, q1,q0, n,ninv ) \
 {                                          \
    umul_ppmm(q1, q0, ninv, p);             \
    add_ssaaaa(q1, q0, q1, q0, p, r);       \
    r -= (q1 + 1) * n;                      \
    if (r >= q0)                            \
     r += n;                                \
    if(r>=n)                                \
     r -= n;                                \
 }
 
static __inline__ 
mp_limb_t NMOD_RED3_pk_func(
  mp_limb_t a_hi,mp_limb_t a_me,mp_limb_t a_lo, 
  const nmod_t mod)
 {
  mp_limb_t t0,t1,r;
  if(a_hi>1)
   {
    umul_ppmm( t0,t1, mod.norm, a_hi );
    a_hi=0;
    add_sssaa( a_hi,a_me,a_lo, t0,t1 );
   }
  // Now a_hi is either 0 or 1
  if( a_hi )
   {
    if( a_me > mod.n )
     a_me -= mod.n;
    a_me -= mod.n;
   }
  // Now result equals 2**64 * a_me + a_lo modulo n
  if( a_me > mod.n )
   a_me -= mod.n;
  #if 0
   r=NMOD_RED2_pk_func( a_me,a_lo, mod.n,mod.ninv );
   return r;
  #else
   NMOD_RED2_pk( a_me,a_lo, t0,t1, mod.n,mod.ninv);
   return a_lo;
  #endif
 }

// result in a_lo. Should be faster than NMOD_RED3
#define NMOD_RED3_pk( a_hi,a_me,a_lo, n,ninv,two_128_mod_n ) \
 {                                                            \
  mp_limb_t RED3_t0,RED3_t1;                                   \
  if(a_hi>1)                                                   \
   {                                                           \
    umul_ppmm( RED3_t0,RED3_t1, two_128_mod_n, a_hi );         \
    a_hi=0;                                                   \
    add_sssaa( a_hi,a_me,a_lo, RED3_t0,RED3_t1 );            \
   }                                                        \
  if( a_hi )                                               \
   {                                                      \
    if( a_me > n )                                       \
     a_me -= n;                                         \
    a_me -= n;                                         \
   }                                                   \
  if( a_me > n )                                       \
   a_me -= n;                                           \
  NMOD_RED2_pk( a_me,a_lo, RED3_t0,RED3_t1, n,ninv);      \
 }

// same as NMOD_RED3_pk, but a_hi<=1
#define NMOD_RED3_pk_easy( a_hi,a_me,a_lo, n,ninv )           \
 {                                                           \
  mp_limb_t t0,t1;                                          \
  if( a_hi )                                               \
   {                                                      \
    if( a_me > n )                                       \
     a_me -= n;                                         \
    a_me -= n;                                         \
   }                                                  \
  if( a_me > n )                                     \
   a_me -= n;                                       \
  NMOD_RED2_pk( a_me,a_lo, t0,t1, n,ninv);        \
 }

#define NMOD_RED3_pk_5arg( a_hi,a_me,a_lo, n,ninv )  \
 {                                                    \
  mp_limb_t t0,t1;                                    \
  NMOD_RED2_pk( a_hi,a_me, t0,t1, n,ninv);            \
  if( a_me > n )                                     \
   a_me -= n;                                       \
  NMOD_RED2_pk( a_me,a_lo, t0,t1, n,ninv);        \
 }

// like count_leading_zeros(), but don't define new variables
#define count_leading_zeros_opt(count, x)                          \
   {                                                                \
    __asm__ ("bsrq %1,%0" : "=r" (count) : "rm" ((mp_limb_t)(x)));  \
    count ^= (mp_limb_t) 63;                                       \
   }

#if 0
z=n_addmod(x,y,n) unrolls into 7 lines without branches:

                         %r12 = n
                         %rbp = x
                         %rbx = y

           	mov    %r12,%rax      rax=n
           	sub    %rbx,%rax      rax=n-y
           	add    %rbp,%rbx      rbx=x+y
           	mov    %rbx,%rsi      rsi=x+y
           	sub    %r12,%rsi      rsi=x+y-n
           	cmp    %rax,%rbp
           	cmovb  %rbx,%rsi      maybe change rsi
#endif

// 2 conditional jumps and two addition-like operators, total 5 statements:
#define n_addmod_pk(t,y,n)                             \
  __asm__                                              \
   (                                                   \
    "addq %1,%q0\n\t"                                  \
    "jc 1f\n\t"                                        \
    "cmpq %2,%0\n\t"                                   \
    "jb 2f\n"                                          \
    "1:subq %2,%0\n"                                   \
    "2:"                                               \
    : "+r" (t)                                         \
    : "rm"( (mp_limb_t)(y) ),"rm" ( (mp_limb_t)(n) )   \
   );

#define n_negmod_opt(t,n) \
 if(t)                    \
  t=n-t;                  \

// r := w*x-y*z modulo n, n >= 2**63
#define WX_MINUS_YZ(r, w,x,y,z, n,ninv )              \
 if(y==0)                                              \
  r=n_mulmod_preinv_4arg(w,x, n,ninv);                 \
 else                                                  \
  {                                                    \
   mp_limb_t _t0,_t1,_t2,_t3=0;                        \
   umul_ppmm(_t1,_t2,n-y,z);                           \
   umul_ppmm( _t0,r, w,x );                            \
   add_sssaa(_t3,_t0,r, _t1,_t2);                       \
   if( _t3 )                                             \
    {                                                     \
     if( _t0 > n )                                       \
      _t0 -= n;                                         \
     _t0 -= n;                                         \
    }                                                 \
   if( _t0 > n )                                     \
    _t0 -= n;                                       \
   NMOD_RED2_pk( _t0,r, _t1,_t2, n,ninv);          \
  }

#endif

#endif
