// This program is part of RAZIN
/*
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009 William Hart
   Copyright 2011 Fredrik Johansson
*/
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL) version 2.1 or

// This file contains a modified fragment of FLINT longlong.h, owned by W.Hart 
//  and F.Johansson

#ifndef LONGLONG__H
#define LONGLONG__H

#if (GMP_LIMB_BITS == 64 && defined (__amd64__))

// slim version of add_sssaaaaaa macro found in FLINT longlong.h
#define add_sssaaa0aa(sh, sm, sl, ah, am, al,    bm, bl)  \
  __asm__ ("addq %8,%q2\n\tadcq %6,%q1\n\tadcq %4,%q0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "N"   (0x0),              \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  \

#endif

#endif
