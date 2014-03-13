#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks different algorithms to find inverse of an integer 
 matrice
'''

import sage.all
import flint_sage as flint
import sys,time

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=30

t_sage=t__inv=t_matz=t_Dixon=0

def random_matrice(dim):
 return sage.all.random_matrix( ZZ, dim, x=-100, y=100 )

def test_with( m ):
 '''
 sage method only called once, do not need to clean cache
 '''
 global t_sage, t__inv, t_matz, t_Dixon
 i=sage.all.identity_matrix( m.nrows() )
 mZ=flint.fmpz_mat( m )
 mQ=flint.fmpq_mat( (Integer(1), mZ) )
 iZ=flint.fmpz_mat( i )
 iQ=flint.fmpq_mat( (Integer(1),iZ) )
 try:
  t0=time.time()
  sS=m.inverse()
 except:
  t_sage += time.time()-t0
  singular_test( mZ, mQ, iZ, iQ )
  return
 t1=time.time()
 t_sage += t1-t0
 Xd,d=mZ.solve_right(iZ)
 t2=time.time()
 t_matz += t2-t1
 Y=mQ.solve_dixon( iQ )
 t3=time.time()
 t_Dixon += t3-t2
 W=mQ.inv()
 t__inv += time.time()-t3

def singular_test( mZ, mQ, iZ, iQ ):
 global t__inv, t_matz, t_Dixon
 t1=time.time()
 drop,zero=mZ.solve_right(iZ)
 t2=time.time()
 X=mQ.solve_dixon( iQ )
 t3=time.time()
 W=mQ.inv()
 t4=time.time()
 t_matz += t2-t1
 t_Dixon += t3-t2
 t__inv += t4-t3

a=sage.all.matrix( 2, [5, 5, 7, 8] )
test_with( a )
a[1,0],a[1,1]=-5,-5
test_with( a )

def print_rez(dim):
 print '%3d   %.2e   %.2e   %.2e   %.2e' % (dim,t_sage,t__inv,t_matz,t_Dixon)

print '        sage       inv        matz       Dixon'
for i in range(33):
 dim=dim0 + 10*i
 t_sage=t__inv=t_matz=t_Dixon=0
 for j in range(9):
  a=random_matrice( dim )
  test_with( a )
 for j in range(dim):
  a[1,j]=0
 test_with( a )
 print_rez(dim)
