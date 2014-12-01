#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks C subroutine fmpz_mat_det_hermitian_decomposition() 
 against fmpz_mat_det_4block() or fmpz_mat_det()

1000-bits entries matrix with varying size taken (5 matrix per one dimension)

20141129 benchmark results in the bottom of the file
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

' 3 subroutines have similar interface, flint.det calls FLINT fmpz_mat_det() '
#count_det=flint.det_20140704
count_det=flint.det
count_hd_det=flint.det_hermitian_decomposition 
        # fmpz_mat_det_hermitian_decomposition

bits=1000
xy_par=1<<bits
profile_on=0

def test_for_dim(n):
 global t0,t1
 t0=t1=0
 if profile_on:
  print 'dim=%s' % n 
 for u in range(5):
  a=random_matrix(ZZ, n, x=-xy_par, y=xy_par)
  test_matrice_mind_time(a)
 print '=',n,t0,t1
 if profile_on:
  print '**************************************************************\n\n'
  
def test_matrice_mind_time(a):
 global t0,t1
 m0=time.time()
 z0=count_det( fmpz_mat(a) )
 m1=time.time()
 t0 += m1-m0
 if profile_on:
  print '--------------------------------------------------------------\n'
 m0=time.time()
 z1=count_hd_det( fmpz_mat(a) )
 m1=time.time()
 if z1 != z0:
  print 'det mismatch: %X != %X' % (z0,z1)
  assert 0
 if profile_on:
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
     FLINT        fmpz_mat_det_hermitian_decomposition()
100 12.8463242054 9.4831635952
200 169.630519152 110.104531765
300 688.357366323 388.271833897
400 1793.7264452 961.077443123
500 3725.16406798 1881.69160128
600 7468.16601181 3322.10169816

Benchmark for Sage 6.1.1-r2 above 
Benchmark for Sage 6.4-r1 below

100 12.0137267113 8.72138190269
200 134.917006969 70.5969481468
300 525.445667267 255.590985775
400 1428.36421037 603.439989328
500 3176.60813093 1186.68403578
600 6105.40285397 2113.61041331
'''
