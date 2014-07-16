// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <flint/flint.h>
#include <flint/longlong.h>
#include "longlong_.h"

#if 0
This program tests macro add_sssaaa0aa

To compile:
gcc add_2.c -march=native -O2 -lm -IC/ulong_extras
#endif

#define size 100
uint64_t A[size];
uint64_t R0[3];
uint64_t R1[3];

void 
random_fill(uint64_t* t,uint64_t siz,uint64_t seed)
 {
  int i;
  double d;
  if( 0==seed ) 
   seed=~seed;
  for(i=siz;i--;)
   {
    d=log2( (double)seed ) + sin( (double)i );
    t[i] = *(uint64_t*)&d;
   }
 }

void slow_way(uint64_t* r,uint64_t* s,uint64_t siz)
 {
  r[0]=r[1]=r[2]=0;
  uint64_t i,j,z=0;
  for(i=0;i<siz;i++)
   {
    j=2*i;
    add_sssaaaaaa( r[2],r[1],r[0], r[2],r[1],r[0], 0, s[j], s[j+1] );
   }
 }

void fast_way(uint64_t* r,uint64_t* s,uint64_t siz)
 {
  r[0]=r[1]=r[2]=0;
  uint64_t i,j;
  for(i=0;i<siz;i++)
   {
    j=2*i;
    add_sssaaa0aa( r[2],r[1],r[0], r[2],r[1],r[0],    s[j], s[j+1] );
   }
 }

int main()
 {
  int i;
  for(i=10;i--;)
   {
    random_fill(A,size,i);
    slow_way(R0,A,size/2);
    fast_way(R1,A,size/2);
    if( (R0[0] != R1[0]) || (R0[1] != R1[1]) || (R0[2] != R1[2]) )
     assert(0);
   }
  printf("Test passed\n");
 }
