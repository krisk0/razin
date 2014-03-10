#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# This program requires flint_sage package
# This program takes time. Does 10 tests. Prints intermediate results after 
#  each of 10 tests. No files created. Run the program and relax

'''
 Benchmark different ways of checking that a matrice is unimodular
 
 Approach 1: count determinant with Sage native subroutine, NTL, PARI/GP,
  LinBox, FLINT. Ups, LinBox only offers unsafe algorithm. LinBox, goodbye
  
 Approach 2: solve for x equation a*x=b where b is a basis vector. a is 
  unimodular iff integer salvation exists for all basis vectors
'''

import sage.all,numpy,time,sys
import flint_sage as flint

sageM=sage.all.matrix
ZZ=sage.all.ZZ

dim=100
I=sage.all.identity_matrix( ZZ, dim )
experiments=10
ways_desc='Sage NTL pari FLINT solve_r'
ways=5

t1=numpy.resize( numpy.array( [], dtype=float ), ways )
tx=t1.copy()
strange_det=0

def random_matrice():
 a=sage.all.random_matrix(ZZ, dim, x=-100, y = 100)
 return a.hermite_form( transformation=True )[1]

def count_det(a,me,i):
 global t1,tx,strange_det
 t0=time.time()
 one=a.determinant(algorithm=me,proof=True)
 t1[i] += time.time()-t0
 if abs(one) != 1:
  print 'bad determinant, method=%s' % me
  print 'a=\n'
 b=sage.all.copy(a)
 b[0,0] += 1
 t0=time.time()
 not_one=b.determinant(algorithm=me,proof=True)
 tx[i] += time.time()-t0
 if abs(not_one)==1:
  strange_det += 1

def det_4_ways(a):
 way=0
 # 'linbox' is the 4th way, disqualified for providing no sure way
 for m in 'padic','ntl','pari': 
  count_det(a,m,way)
  a._clear_cache()
  way += 1

def solve_right(a,way):
 global t1,tx,strange_det
 t0=time.time()
 for x in range(dim):
  try:
   one=a.solve_right( I[x] ).denominator()
  except:
   one=0
  if abs(one) != 1:
   t1[way] += time.time()-t0
   print 'bad determinant, solve_right' % me
   print 'a=\n',a
   break
 if abs(one)==1:
  t1[way] += time.time()-t0
 b=sage.all.copy(a)
 b[0,0] += 1
 t0=time.time()
 for x in range(dim):
  try:
   not_one=b.solve_right( I[x] ).denominator()
  except:
   not_one=0
  if not_one != 1:
   tx[way] += time.time()-t0
   break
 if abs(not_one)==1:
  tx[way] += time.time()-t0
  strange_det += 1

def det_flint(a,way):
 global t1,tx,strange_det
 aF=flint.fmpz_mat( a )
 t0=time.time()
 one=flint.det(aF)
 t1[way] += time.time()-t0
 if abs(one) != 1:
  print 'bad determinant, flint'
  print 'a=\n',a
 b=sage.all.copy(a)
 b[0,0] += 1
 bF=flint.fmpz_mat( b )
 t0=time.time()
 not_one=flint.det(bF)
 tx[way] += time.time()-t0
 if abs(not_one)==1:
  strange_det += 1

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
 if strange_det:
  # strange_det should be multiple of ways
  print 'Warped matrice is also unimodular: %s times' % strange_det

for x in range(experiments):
 a=random_matrice()
 det_4_ways(a)
 det_flint(a,3)
 solve_right(a,4)
 for x in range(ways):
  sys.stdout.write( '%.4e ' % t1[x] )
 print

pretty_print_result()
