// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/
// Copyright Денис Крыськов 2014

// This file contains modified FLINT function n_pow()

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

mp_limb_t n_pow_speedup(mp_limb_t n, ulong exp)
// This function does one multiplication less than original
{
   ulong i;
   mp_limb_t res;

   res = n;
   for (i = exp-1; i--; )
      res *= n;

   return res;
}
