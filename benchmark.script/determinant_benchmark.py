#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# This program requires flint_sage package
# This program takes time. Does 10 tests. Prints intermediate results after 
#  each of 10 tests. No files created. Run the program and relax

'''
 Benchmark determinant calculation

 Result for random U (not unimodular) at the end of file

 pari terribly slow, flint champion, sage near
'''

import sage.all,numpy,time,sys,os
import flint_sage as flint

sageM=sage.all.matrix
ZZ=sage.all.ZZ

dim=100
I=sage.all.identity_matrix( ZZ, dim )
experiments=3
ways_desc='Sage NTL pari FLINT'
ways=4
loud=os.environ.get('loud')

t1=numpy.resize( numpy.array( [], dtype=float ), ways )
tx=t1.copy()

def random_matrice():
 if 0:
  a=sage.all.random_matrix(ZZ, dim, x=-1<<128, y = 1<<128)
  return a.hermite_form( transformation=True )[1]
 else:
  return sage.all.random_matrix(ZZ, dim, x=-1<<128, y = 1<<128)

def count_det(a,me,i):
 global t1,tx
 t0=time.time()
 one=a.determinant(algorithm=me,proof=True)
 t1[i] += time.time()-t0
 if loud:
  print one
 b=sage.all.copy(a)
 b[0,0] += 1
 t0=time.time()
 not_one=b.determinant(algorithm=me,proof=True)
 tx[i] += time.time()-t0
 if loud:
  print not_one

def det_4_ways(a):
 way=0
 # 'linbox' is the 4th way, disqualified for providing no sure way
 for m in 'padic','ntl','pari': 
  count_det(a,m,way)
  a._clear_cache()
  way += 1

def det_flint(a,way):
 global t1,tx
 aF=flint.fmpz_mat( a )
 t0=time.time()
 one=flint.det(aF)
 t1[way] += time.time()-t0
 if loud:
  print one
 b=sage.all.copy(a)
 b[0,0] += 1
 bF=flint.fmpz_mat( b )
 t0=time.time()
 not_one=flint.det(bF)
 tx[way] += time.time()-t0
 if loud:
  print not_one

def chew_desc( d ):
 r,l=d.split(' '),0
 for x in r:
  s=len(x)
  if s>l:
   l=s
 return r,l

def pretty_print_result():
 print
 desc,max_l=chew_desc( ways_desc )
 print (' '*max_l)+'     U      warped U'
 print ('-'*max_l)+' --------- ---------'
 fmt='%'+str(max_l)+'s  %.2e  %.2e'
 for i in range(ways):
  print fmt % (desc[i],t1[i],tx[i])

for x in range(experiments):
 a=random_matrice()
 det_4_ways(a)
 det_flint(a,3)
 #solve_right(a,4)
 for x in range(ways):
  sys.stdout.write( '%.4e ' % t1[x] )
 print

pretty_print_result()
'''
            U      warped U
------- --------- ---------
   Sage  7.77e-01  7.86e-01
    NTL  1.14e+00  1.13e+00
   pari  1.89e+03  1.46e+03
  FLINT  7.31e-01  7.30e-01
'''
