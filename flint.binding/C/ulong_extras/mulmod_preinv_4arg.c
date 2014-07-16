// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2009, 2010, 2012 William Hart

******************************************************************************/
// Copyright Денис Крыськов 2014

// This file contains a piece of FLINT mulmod_preinv.c owned by W.Hart

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

mp_limb_t
n_mulmod_preinv_4arg(mp_limb_t a, mp_limb_t b, mp_limb_t n, mp_limb_t ninv)
//castrated FLINT subroutine n_mulmod_preinv(): norm parameter removed
{
    mp_limb_t q0, q1, r, p;

    /* multiply */
    umul_ppmm(p, r, a, b);
    
    /* reduce mod n */
    umul_ppmm(q1, q0, ninv, p);
    add_ssaaaa(q1, q0, q1, q0, p, r);

    r -= (q1 + 1) * n;

    if (r >= q0)
     r += n;

    return (r < n ? r : r - n);
}
