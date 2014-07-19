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

#include <flint/nmod_vec.h>

#if (GMP_LIMB_BITS == 64 && defined (__amd64__))

// slim version of add_sssaaaaaa macro found in FLINT longlong.h
#define add_sssaaa0aa(sh, sm, sl, ah, am, al,    bm, bl)  \
  __asm__ ("addq %8,%q2\n\tadcq %6,%q1\n\tadcq %4,%q0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "N"   (0x0),              \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  \

/*
#define add_sss0aa0aa(sh, sm, sl,     am, al,    bm, bl)  \
  __asm__ ("addq %6,%q2\n\tadcq %4,%q1\n\tadcq $0x0,%q0"     \
       : "=r" (sh), "=&r" (sm), "=&r" (sl)                  \
       : "0"  (( ???                                                 \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  \
*/

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

// result in r
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
  const nmod_t mod, 
  mp_limb_t two_128_mod_n)
 {
  mp_limb_t t0,t1,r;
  if(a_hi>1)
   {
    umul_ppmm( t0,t1, two_128_mod_n, a_hi );
    a_hi=0;
    add_sssaaa0aa( a_hi,a_me,a_lo, a_hi,a_me,a_lo, t0,t1 );
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
  mp_limb_t t0,t1;                                             \
  if(a_hi>1)                                                   \
   {                                                           \
    umul_ppmm( t0,t1, two_128_mod_n, a_hi );                   \
    a_hi=0;                                                   \
    add_sssaaa0aa( a_hi,a_me,a_lo, a_hi,a_me,a_lo, t0,t1 );  \
   }                                                        \
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

  40063a:	4c 89 e0             	mov    %r12,%rax      rax=n
  40063d:	48 29 d8             	sub    %rbx,%rax      rax=n-y
  400640:	48 01 eb             	add    %rbp,%rbx      rbx=x+y
  400643:	48 89 de             	mov    %rbx,%rsi      rsi=x+y
  400646:	4c 29 e6             	sub    %r12,%rsi      rsi=x+y-n
  400649:	48 39 c5             	cmp    %rax,%rbp
  40064c:	48 0f 42 f3          	cmovb  %rbx,%rsi      if condition holds, rsi=x+y
                                                      else leave rsi be x+y-n
#endif

// 2 conditional jumps and two addition-like operators:
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

#endif

#endif
