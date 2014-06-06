// This program is part of RAZIN
/******************************************************************************
   
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart

******************************************************************************/
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "flint/flint.h"

static __inline__ int
n_is_strong_probabprime2_preinv_speedup(mp_limb_t n, mp_limb_t ninv, mp_limb_t a,
                                mp_limb_t d)
/*
this subroutine does Miller-Rabin test and returns positive iff test passes

hacked by Денис Крыськов to count n-1 once
*/
 {
    mp_limb_t t = d;
    mp_limb_t y;

    y = n_powmod2_ui_preinv(a, t, n, ninv);

    if (y == UWORD(1))
        return 1;
    t <<= 1;

    d = n-1;        // Денис Крыськов was here
    while ((t != d) && (y != d))  // and here
     {
        y = n_mulmod2_preinv(y, y, n, ninv);
        t <<= 1;
     }
    return y == d;            // and here
 }
