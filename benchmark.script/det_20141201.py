#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks C subroutine fmpz_mat_det_hermitian_decomposition() 
 against 2 friends

1000-bits entries matrix with varying size taken (3 to 5 matrix per one 
 dimension)
'''

import sage.all
import flint_sage as flint
import sys,time,numpy
from test_nmod_HNF import unimodular,unimodular_triU

matrix=sage.all.matrix
ZZ=sage.all.ZZ
write=sys.stdout.write

randint=sage.all.randint
random_matrix=sage.all.random_matrix
fmpz_mat=flint.fmpz_mat

algo_desc=['NTL','FLINT','RAZIN']
bits=1000

xy_par=1<<bits
algo_count=len(algo_desc)
time_array=func_array=None

def test_matrice(a):
 '''
  Run algorithms in a random order
  
  Abort if results differ
  
  Store time in time_array
 '''
 global time_array
 p=random_permutation(algo_count)
 prev_det=None
 for i in range(algo_count):
  t0=time.time()
  curr_det=func_array[ p[i] ](a)
  t1=time.time()
  if None==prev_det:
   prev_det=curr_det
  else:
   if prev_det!=curr_det:
    print "Algorithms %s and %s gave different results on matrice" \
     % (algo_desc[ p[i-1] ], algo_desc[ p[i] ])
    print a
    sys.exit(1)
  time_array[ p[i] ] += t1-t0

def det_NTL(a):
 return a.determinant(algorithm='ntl',proof=True)
 '''
 No need to clear cache
 
 r=a.determinant(algorithm='ntl',proof=True)
 a._clear_cache()
 return r
 '''

def det_FLINT(a):
 return flint.det( fmpz_mat(a) )

def det_RAZIN(a):
 return flint.det_hermitian_decomposition( fmpz_mat(a) )

def random_permutation(n):
 r=numpy.resize( numpy.array([],dtype=int), algo_count )
 for i in range(n):
  r[i]=i
 for i in range(n-1):
  j=randint( i+1, n-1 )
  r[i],r[j]=r[j],r[i]
 return r

def print_time(n, i):
 write('%3s ' % n)
 for k in range(algo_count):
  write( pretty_float(time_array[k]/i) )
 print

def pretty_float(x):
 r='%8f' % x
 if len(r)>8:
  p=r.find('.')
  r=('%'+str(p)+'.'+str(7-p)+'f') % x
 return r+' '

def test_for_dim(n,i):
 global time_array
 for k in range(algo_count):
  time_array[k]=0
 for k in range(i):
  a=random_matrix(ZZ, n, x=-xy_par, y=xy_par)
  test_matrice(a)
 print_time(n,i)
  
def chew_desc(d, p):
 global time_array,func_array
 fmt='%'+('%s' % p)+'s'
 write('  ')
 for i in range(algo_count):
  write( fmt % algo_desc[i] )
 print
 time_array=numpy.resize( numpy.array([]), algo_count )
 func_array=numpy.array([det_NTL,det_FLINT,det_RAZIN],dtype=object)

sage.all.set_random_seed('20141201')
chew_desc(algo_desc,9)

for i in range(10):
 test_for_dim(10+i*10,5)
for i in range(5):
 test_for_dim(200+i*100,3)

print '\ntest passed'

'''
Intel(R) Core(TM) i5-2500K CPU @ 3.30GHz
        NTL    FLINT    RAZIN
 10 0.010621 0.001683 0.007121 
 20 0.020692 0.039772 0.028679 
 30 0.066989 0.067361 0.071621 
 40 0.163061 0.146227 0.143874 
 50 0.333710 0.277872 0.253775 
 60 0.606659 0.491468 0.401507 
 70 1.013421 0.803321 0.599855 
 80 1.596735 1.205350 0.893154 
 90 2.396125 1.733157 1.213826 
100 3.459566 2.423262 1.638669 
200 42.74575 27.64102 13.49951 
300 193.8454 107.0062 49.17864 
400 577.7632 289.3752 114.4222 
500 1359.483 643.6196 226.2774 
600 2753.817 1242.382 400.9209 
'''
