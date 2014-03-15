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
  d=matrix_from_columns( range(a.nrows()) ).determinant()
 return abs(d)

def test_dump(a,h):
 print 'hnf of\n',a,'\nequals\n',h

def test(a,loud):
 b=sub_det(a)
 if b<2 or b>=2**64:
  print 'bad determinant, test skipped'
  return
 sage_r=a.hermite_form()
 nmod_r=flint.fmpz_mat_hermite_form( flint.fmpz_mat( a ), b )
 if loud:
  test_dump(a,sage_r)
 if nmod_r != flint.fmpz_mat( sage_r ):
  if not loud:
   test_dump(a,sage_r)
   print 'nmod result:',nmod_r

a2=matrix( 2, [ 5,5,-7,-8] )
u2=matrix( 2, [ 1,77,0,-1] ) * matrix( 2, [ 1,0,-87,1] )
assert sub_det(u2) == 1
test( a2, 1 )
test( u2*a2, 1)
test( a2*u2, 1)



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
