#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys

dim=100
ZZ=sage.all.ZZ

def random_matrice():
 a=sage.all.random_matrix(ZZ, dim, x=-100, y = 100)
 return a

def test_with(m, precook):
 i=sage.all.identity_matrix( m.nrows() )
 mF=flint.fmpz_mat( m )
 # leave i as-is or convert to fmpz_mat
 if precook:
  iF=flint.fmpz_mat( i )
 else:
  iF=i
 try:
  s=m.solve_right(i)
 except:
  Xd,d=mF.solve_right(iF)
  assert d==0
  sys.write('-')
  return
 # matrix non-singular, can do arithmetic
 Xd,d=mF.solve_right(iF)
 # ups, where is that Matrix_integer_dense?
 #Sd=sage.all.Matrix_integer_dense( d*s )
 Sd=sage.all.matrix( ZZ, d*s )
 Sd=flint.fmpz_mat( Sd )
 if Xd != Sd:
  print '\n test failed,\n m=\n%s' % m
  print 's=%s' % s
  print 'd=%d' % d
  print 'd=%Xd' % Xd
  sys.exit(1)
 #sys.stdout.write('+')
 sys.stdout.write( ' %s' % d )

a=sage.all.matrix( 2, [5, 5, 7, 8] )
test_with( a, 0 ),test_with( a, 1 )

sys.stdout.write('\n+')
for i in range(10):
 a=random_matrice()
 # test both forms of call
 test_with( a, 0 ),test_with( a, 1 )
 print

print '\ntest passed'
