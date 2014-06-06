// This program is part of RAZIN

// Division-free 1-limb GCD code inspired by A. Baid version of ulong_extras/gcd.c
/******************************************************************************
   
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "flint/flint.h"

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
 
mp_limb_t
n_gcd_odd_odd(mp_limb_t x,mp_limb_t y)
/*
x,y odd
*/
 {
  register mp_limb_t s;
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
