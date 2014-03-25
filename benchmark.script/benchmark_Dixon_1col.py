#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks Sage .solve_right() and flint_sage .solve_right() and 
 .solve_Dixon()
 
This program outputs extra lines such as ++++-
 It is a feature and not bug
'''

import sage.all
import flint_sage as flint
import sys,time

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=20
t_sage=t_matz=t_matq=0

def random_matrice(dim):
 return sage.all.random_matrix( ZZ, dim, x=-100, y=100 )

def test_with( m ):
 global t_sage, t_matz, t_matq
 i=sage.all.identity_matrix( m.nrows() ).column(1).column()
 mZ=flint.fmpz_mat( m )
 mQ=flint.fmpq_mat( (Integer(1), mZ) )
 iZ=flint.fmpz_mat( i )
 iQ=flint.fmpq_mat( (Integer(1),iZ) )
 singular=0
 t0=time.time()
 try:
  sS=m.solve_right(i)
 except:
  write('-')
  singular=1
 t1=time.time()
 Xd,d=mZ.solve_right(iZ)
 t2=time.time()
 Y=mQ.solve_dixon( iQ )
 t3=time.time()
 t_sage += t1-t0
 t_matz += t2-t1
 t_matq += t3-t2
 if not singular:
  write('+')

def print_rez(dim):
 print '\n%3d   %.2e    %.2e     %.2e' % (dim,t_sage,t_matz,t_matq)
 sys.stdout.flush()

print '        sage    solve_right      Dixon'
for i in range(33):
 t_sage=t_matz=t_matq=0
 dim=dim0 + 10*i
 for j in range(77):
  a=random_matrice( dim )
  test_with( a )
 for j in range(dim):
  a[1,j]=0
 test_with( a )
 print_rez(dim)
