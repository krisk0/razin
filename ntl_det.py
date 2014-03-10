#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

'''
 Benchmark NTL detereminant with proof=True and proof=False
'''

import sage.all,numpy,time,sys

sageM=sage.all.matrix
ZZ=sage.all.ZZ

dim=100
experiments=15
ways_desc='safe unsafe'
ways=2

t1=numpy.resize( numpy.array( [], dtype=float ), ways )
tx=t1.copy()
strange_det=unsafe_fail=0

def random_matrice():
 a=sage.all.random_matrix(ZZ, dim, x=-100, y = 100)
 u=a.hermite_form( transformation=True )[1]
 return u._ntl_()

def count_det(a,me,i):
 global t1,tx,strange_det
 t0=time.time()
 one=int(a.determinant(deterministic=me))
 end_t=time.time()
 print 'method=%s a start time %f end time %f' % (i,t0,end_t)
 t1[i] += end_t-t0
 if abs(one) != 1:
  print 'bad determinant, method=%s' % me
  print 'a=\n'
 b=sage.all.copy(a)
 b[0,0] += 1
 t0=time.time()
 not_one=int(b.determinant(deterministic=me))
 end_t=time.time()
 print 'method=%s b start time %f end time %f' % (i,t0,end_t)
 tx[i] += end_t-t0
 if abs(not_one)==1:
  strange_det += 1
 return (one,not_one)

def det_2_ways(a):
 global unsafe_fail
 safe_r=count_det(a,True,0)
 unsafe_r=count_det(a,False,1)
 if safe_r != unsafe_r:
  unsafe_fail += 1

def chew_desc( d ):
 r,l=d.split(' '),0
 for x in r:
  s=len(x)
  if s>l:
   l=s
 return r,l

def pretty_print_result():
 print 't1=',t1
 print 'tx=',tx
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
 if unsafe_fail:
  print 'unsafe method failed:',unsafe_fail

for x in range(experiments):
 a=random_matrice()
 det_2_ways(a)
 for x in range(ways):
  sys.stdout.write( '%.4e ' % t1[x] )
 print

pretty_print_result()
