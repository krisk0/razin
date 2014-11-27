#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks C subroutine fmpz_mat_det_hermitian_decomposition() 
 against fmpz_mat_det_4block() or fmpz_mat_det()

1000-bits entries matrix with varying size taken (5 matrix per one dimension)

some results are in the bottom of the file
'''

import sage.all
import flint_sage as flint
import sys,time
from test_nmod_HNF import unimodular,unimodular_triU

Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

randint=sage.all.randint
random_matrix,identity_matrix=sage.all.random_matrix,sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat

#count_det=flint.det_20140704
count_det=flint.det
count_hd_det=flint.det_hermitian_decomposition

bits=1000
xy_par=1<<bits

def test_for_dim(n):
 global t0,t1
 t0=t1=0
 print 'dim=%s' % n 
 for u in range(5):
  a=random_matrix(ZZ, n, x=-xy_par, y=xy_par)
  test_matrice_mind_time(a)
 print '=',n,t0,t1
 print '**************************************************************\n\n'
  
def test_matrice_mind_time(a):
 global t0,t1
 m0=time.time()
 z0=count_det( fmpz_mat(a) )
 m1=time.time()
 t0 += m1-m0
 print '--------------------------------------------------------------\n'
 m0=time.time()
 z1=count_hd_det( fmpz_mat(a) )
 m1=time.time()
 if z1 != z0:
  print 'det mismatch: %X != %X' % (z0,z1)
  assert 0
 print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
 t1 += m1-m0

sage.all.set_random_seed('20141126')

t0=t1=0
for i in range(10):
 test_for_dim(10+i*10)
for i in range(5):
 test_for_dim(200+i*100)

print '\ntest passed'

'''
My method appears to be asymptotically better

100 12.8778271675 9.85142207146
200 172.378462076 110.319196224
300 689.339381933 404.530629873
400 1796.41489267 970.717617989
500 3727.82309914 1934.79675388
600 7395.46039391 3358.61122966
'''
