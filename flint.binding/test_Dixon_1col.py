#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program solves a system of linear equations in different ways, thus testing 
 wrappers for fmpz_mat_solve, fmpq_mat_solve_dixon

if watch_canonical=1, then print some information if raw_str() returns different
 information (which can only happen when one of the matrix is not in canonical 
 form)
'''

import sage.all
import flint_sage as flint
import sys

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=15
watch_canonical=1

def random_matrice(dim):
 return sage.all.random_matrix( ZZ, dim, x=-100, y=100 )

def test_with( m ):
 i=sage.all.identity_matrix( m.nrows() ).column(1).column()
 assert i.nrows()==m.nrows()
 mZ=flint.fmpz_mat( m )
 mQ=flint.fmpq_mat( (Integer(1), mZ) )
 iZ=flint.fmpz_mat( i )
 iQ=flint.fmpq_mat( (Integer(1),iZ) )
 try:
  sS=m.solve_right(i)
 except:
  write('-')
  singular_test( mZ, mQ, iZ, iQ )
  return
 Xd,d=mZ.solve_right(iZ)
 S_by_d_sage=sage.all.matrix( ZZ, d*sS )
 S_by_d=flint.fmpz_mat( S_by_d_sage )
 if Xd != S_by_d:
  print 'integer solve_right() failed, m=\n',m
  sys.exit(1)
 X=flint.fmpq_mat( (d, Xd) )
 Y=mQ.solve_dixon( iQ )
 assert X == Y
 write('+')
 if watch_canonical:
  ' check if raw_str() of X and Y are equal '
  Xr,Yr=X.raw_str(),Y.raw_str()
  if Xr != Yr:
   print 'm=\n',m
   print 'Xr=%s' % Xr
   print 'Yr=%s' % Yr

def singular_test( mZ, mQ, iZ, iQ ):
 drop,zero=mZ.solve_right(iZ)
 assert zero==0
 X=mQ.solve_dixon( iQ )
 assert X==None
 W=mQ.inverse()
 assert W==None

a=sage.all.matrix( 2, [5, 5, 7, 8] )
test_with( a )
a[1,0],a[1,1]=-5,-5
test_with( a )

sys.stdout.write('\n')
for i in range(10):
 dim=dim0 + 5*i
 for j in range(7):
  a=random_matrice( dim )
  test_with( a )
 for j in range(dim):
  a[1,j]=0
 test_with( a )
 print
 
print '\ntest passed'
