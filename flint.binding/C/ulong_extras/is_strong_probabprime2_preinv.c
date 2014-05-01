/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************
   
    Copyright (C) 2008, Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Денис Крыськов

     Borrowed some code from Abhinav Baid re-implemented n_gcd()

******************************************************************************/

#include "flint/flint.h"
//#include "ulong_extras.h"

#if (defined (__amd64__) || defined (__i386__) || defined (__i486__)) 
mp_limb_t
n_gcd_odd_even(mp_limb_t x,mp_limb_t y)
/*
x is odd, y is possibly even
 
Inspired by division-free n_gcd() subroutine by Abhinav Baid
*/
 {
  register mp_limb_t s;
  count_trailing_zeros(s, y);
  y >> s;
  while(x!=y)
   {
    if(x<y)
     {
       y-=x;
       count_trailing_zeros(s, y);
       y >>= s;
     }
    else
     {
       x-=y;
       count_trailing_zeros(s, x);
       x >>= s;
     }
   }
  return x;
 }
#endif

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
