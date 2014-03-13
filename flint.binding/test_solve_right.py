#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
flint solve_right() on matrice with zero row crashes... no longer
'''

import sage.all
import flint_sage as flint
import sys

dim=10
ZZ=sage.all.ZZ
write=sys.stdout.write
experience_flint_crash=1      # if 1 then attempt to crash by using singular matrice

def random_matrice():
 return sage.all.random_matrix( ZZ, dim, x=-100, y=100 )

def test_with(m, precook):
 i=sage.all.identity_matrix( m.nrows() )
 mF=flint.fmpz_mat( m )
 # leave i as-is or convert to fmpz_mat
 if precook:
  iF=flint.fmpz_mat( i )
 else:
  iF=i
 try:
  s=m.solve_right(i)    # yep, this line works twice for every matrice
 except:
  print '\n singular case, dim=%s, calling flint solve_right()' % m.nrows()
  try:
   Xd,d=mF.solve_right(iF) 
  except:
   print 'singular matrice, exception in flint solve_right()'
   return
  print '\n singular case, no exception in flint solve_right()'
  assert d==0
  write('0 ')
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
 if 0:
  write('+')
 else:
  write( ' %s' % d )

a=sage.all.matrix( 2, [5, 5, 7, 8] )
test_with( a, 0 ),test_with( a, 1 )
if experience_flint_crash:
 a[1,0],a[1,1]=-5,-5
 print 'Singular matrice:\n',a
 test_with( a, 0 ),test_with( a, 1 )

sys.stdout.write('\n+')
for i in range(10):
 a=random_matrice()
 # test both forms of call
 test_with( a, 0 ),test_with( a, 1 )
 print
 if experience_flint_crash and i==0:
  ' test matrice with zero row'
  for j in range(dim):
   a[0,j]=0
  test_with( a, 0 ),test_with( a, 1 )
 if experience_flint_crash and i==1:
  ' test matrice with zero column'
  for j in range(dim):
   a[j,0]=0
  test_with( a, 0 ),test_with( a, 1 )
 
print '\ntest passed'
