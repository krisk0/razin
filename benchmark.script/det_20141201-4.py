#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks C subroutine fmpz_mat_det_hermitian_decomposition() 
 against 3 friends

Sage benchmark is only testing speed of Python code calling FLINT fmpz_mat_det()
Therefore you better run det_20141201.py instead of this program

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

algo_desc=['Sage','NTL','FLINT','RAZIN']
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
  #print 'Algorithm [%s]->%s spent %s' % (p[i],algo_desc[ p[i] ],t1-t0)
  time_array[ p[i] ] += t1-t0

def det_Sage(a):
 r=a.determinant(proof=True)
 a._clear_cache()
 return r

def det_NTL(a):
 r=a.determinant(algorithm='ntl',proof=True)
 a._clear_cache()
 return r

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
 #print 'random permutation:',list(r)
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
 func_array=numpy.array([det_Sage,det_NTL,det_FLINT,det_RAZIN],dtype=object)

sage.all.set_random_seed('20141201')
chew_desc(algo_desc,9)

for i in range(10):
 test_for_dim(10+i*10,5)
for i in range(5):
 test_for_dim(200+i*100,3)

print '\ntest passed'

'''
Sage 6.7-r1 (sage-on-gentoo 23 Jun 2015)
-O2 -march=native 
Haswell 
Intel(R) Core(TM) i5-4460  CPU @ 3.20GHz
cache size      : 6144 KB

       Sage      NTL    FLINT    RAZIN
 10 0.001307 0.003391 0.001321 0.005996 
 20 0.031016 0.018190 0.031142 0.023646 
 30 0.054067 0.057901 0.054063 0.057546 
 40 0.122729 0.138530 0.121958 0.116879 
 50 0.233780 0.282431 0.230866 0.212972 
 60 0.410436 0.505625 0.407519 0.338244 
 70 0.661685 0.849899 0.655197 0.517989 
 80 0.979294 1.324821 0.975092 0.768370 
 90 1.425628 1.984981 1.427683 1.054405 
100 1.978595 2.835268 1.965328 1.420375 
200 22.39407 33.92717 22.52225 11.24674 
300 89.56995 154.2797 86.28581 59.01810 
400 236.1032 455.5518 229.6776 141.3112 
500 519.0208 1067.974 514.2779 271.8620 
'''
