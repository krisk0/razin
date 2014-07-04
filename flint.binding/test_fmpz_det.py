#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests subroutine fmpz_mat_det_Kryskov()
'''

import sage.all
import flint_sage as flint,sys

randint=sage.all.randint
identity_matrix=sage.all.identity_matrix
ZZ=sage.all.ZZ

sage.all.set_random_seed('20140704')

def max_degree(n):
 k,m=1,n
 assert n>1
 while 1:
  m *= n
  if m>2**64:
   return k
  k += 1

def test_for_dim( dim ):
 for i in range(100):
  m=sage.all.random_matrix( sage.all.ZZ, dim, x=0, y=100 )
  test_m(m)
  
def test_m(m):
  a=flint.fmpz_mat( m )
  dG = flint.det(a)
  dB = flint.det_20140704( a )
  if dG==dB:
   sys.stdout.write('.')
   return
  print 'test failed, det good/bad=%s/%s' % (dG,dB)
  print m
  assert 0

x=identity_matrix(ZZ, 2)
test_m(x)
x[1,1]=-2
test_m(x)
x=identity_matrix(ZZ, 4)
test_m(x)
x[0,2]=20140704
test_m(x)
x[1,3]=80
test_m(x)

for d in range(5,15):
 sys.stdout.write('%X ' % d)
 test_for_dim( d )

print '\ntest passed'
