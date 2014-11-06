#!/usr/bin/python2 -B -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys,numpy

'''
Test tmod_mat_invert_transpose()

For random matrice A with odd determinant, check that A.T*invert_transpose(A) 
 equals identity matrice modulo 2**64

Todo: test degenerate case, too
'''

ZZ=sage.all.ZZ
identity_matrix=sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat
matrix=sage.all.matrix
write=sys.stdout.write
randint=sage.all.randint
invert_transpose=flint.tmod_mat_invrt_trnspse

from test_nmod_HNF import unimodular,unimodular_triL,unimodular_triU,\
 unimodular_diag

w_modulo=2**64

def do_with(A):
 B=invert_transpose( fmpz_mat(A) )
 dim=A.nrows()
 if identity_matrix(dim) != (B.T*A % 2**64):
  print 'error for A=\n',A
  print 'B=\n',B
  assert 0
 write('+')

def odd_hermite(dim,m):
 a=unimodular_triU(dim,m)
 for i in range(dim):
  a[i,i]=1+2*randint(1,m)
 return a

def check_dim(n):
 a=identity_matrix(n)
 do_with(a)
 a[0,0]=1+2*randint(-99,99)
 do_with(a)
 a[n-1,n-1]=1+2*randint(-99,99)
 do_with(a)
 a[0,0],  a[0,n-1]  =0,1
 a[n-1,0],a[n-1,n-1]=1,0
 do_with(a)
 for i in range(5):
  do_with( unimodular_triL(n,9) )
  do_with( unimodular_triU(n,9) )
  do_with( odd_hermite(n,4) )
  do_with( odd_hermite(n,4).T )
  do_with( unimodular(n)*odd_hermite(n,99) )
  do_with( unimodular(n)*unimodular_diag(n)*odd_hermite(n,99) )

sage.all.set_random_seed('20141102')

do_with( matrix(ZZ,2,[1,0,2,1]) )
do_with( matrix(ZZ,2,[1,0,99,1]) )
do_with( matrix(ZZ,2,[1,99,0,5]) )
do_with( matrix(ZZ,2,[1,0,99,5]) )
do_with( matrix(ZZ,2,[3,0,99,1]) )
do_with( matrix(ZZ,2,[3,0,99,5]) )

for dim in range(2,34):
 check_dim(dim)
 write('.')

print '\ntest passed'
