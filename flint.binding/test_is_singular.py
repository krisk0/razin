#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests subroutine fmpz_mat_is_singular()
'''

import sage.all
import flint_sage as flint,sys

is_singular=flint.fmpz_mat_is_singular_wr

sage.all.set_random_seed('20140807')

def test_for_dim( dim ):
 for i in range(100):
  m=sage.all.random_matrix( sage.all.ZZ, dim, x=0, y=100 )
  test_m(m,0)
  make_singular(m,dim)
  test_m(m,1)

def make_singular(m,n):
 if n==2:
  m[1,0]=7*m[0,0]
  m[1,1]=7*m[0,1]
 else:
  l=n-1;
  for i in range(n):
   m[l,i]=m[0,i]-m[1,i]
  
def test_m(m,s):
  a=flint.fmpz_mat( m )
  i = is_singular(a)
  if s and i:
   return
  d = flint.det_20140704( a )
  if d:
   assert i==0
  else:
   assert i

for d in range(2,15):
 sys.stdout.write('%X ' % d)
 test_for_dim( d )

print '\ntest passed'
