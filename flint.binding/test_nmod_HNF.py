#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests fmpz_mat_hermite_form() function, which is a wrapper above
  nmod_mat.hermite_form()
'''

import sage.all
import flint_sage as flint
import random,sys
Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

def sub_det(a):
 if a.nrows() == a.ncols():
  d=a.determinant()
 else:
  d=a.matrix_from_columns( range(a.nrows()) ).determinant()
 return abs(d)

def test_dump(a,h):
 print 'hnf of\n',a,'\nequals\n',h

def test(a,loud):
 b=sub_det(a)
 if b<2 or b>=2**63:
  print 'bad determinant, test skipped'
  return
 b *= 2
 sage_r=a.hermite_form()
 nmod_r=flint.fmpz_mat_hermite_form( flint.fmpz_mat( a ), b )
 if loud:
  test_dump(a,sage_r)
  print 'nmod result:',nmod_r
 if nmod_r != flint.fmpz_mat( sage_r ):
  if not loud:
   test_dump(a,sage_r)
   print 'nmod result:',nmod_r
  print 'test failed'
  sys.exit(1)

def m(x):
 dim=int( len(x)**.5 )
 return matrix( dim, x )

a2=m( [ 5,5,-7,-8] )

h=flint.fmpz_mat_hermite_form( flint.fmpz_mat( a2 ), Integer(13) )
print 'h=',h
h=flint.fmpz_mat_hermite_form( flint.fmpz_mat( a2 ), Integer(10) )
print 'h=',h
assert 0

u2=m( [ 1,77,0,-1] ) * m( [ 1,0,-87,1] )
assert sub_det(u2) == 1
test( a2, 1 )
test( u2*a2, 1)
test( a2*u2, 1)

v2=m( [ 1,67,0,-1] )
b2=a2 * u2 * v2 * m( [ 3,0,30,36] )

test( b2, 1 )
test( b2 * m( [-1, 17, 0, 1 ] ), 1 )

x=y=2**64-1
z=2**64

print 'test passed'
sys.exit(0)

def test_serie(dim,vol,loud):
 for i in range(vol):
  test ( sage.all.random_matrix(ZZ, dim, x=x, y=y), random.randint(2,z), loud )

test_serie(3,3,1)
test_serie(33,333,0)

print 'test passed'
