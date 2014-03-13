#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=4

def random_matrice(dim):
 return sage.all.random_matrix( ZZ, dim, x=-100, y=100 )

def sage_ZZ_to_flint_QQ( d, a ):
 return flint.fmpq_mat( (Integer(d), flint.fmpz_mat( a ) ) )

def test_with( m ):
 i=sage.all.identity_matrix( m.nrows() )
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
 W=mQ.inv()
 assert X == Y
 assert Y == W
 write('+')

def singular_test( mZ, mQ, iZ, iQ ):
 drop,zero=mZ.solve_right(iZ)
 assert zero==0
 X=mQ.solve_dixon( iQ )
 assert X==None
 W=mQ.inv()
 assert W==None

a=sage.all.matrix( 2, [5, 5, 7, 8] )
test_with( a )
a[1,0],a[1,1]=-5,-5
print 'Singular matrice:\n',a
test_with( a )

sys.stdout.write('\n+')
for i in range(10):
 dim=dim0 + 2*i
 a=random_matrice( dim )
 test_with( a )
 for j in range(dim):
  a[1,j]=0
 test_with( a )
 print
 
print '\ntest passed'
