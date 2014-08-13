#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program benchmarks my subroutine fmpz_mat_det_4block() against FLINT 
 fmpz_mat_det()
'''

import sage.all
import flint_sage as flint,sys,time

randint=sage.all.randint
identity_matrix=sage.all.identity_matrix
ZZ=sage.all.ZZ

sage.all.set_random_seed('20140704')
t_FLINT=t_RAZIN=None

def print_statistic(d):
 if t_RAZIN:
  print '%3d: %5f %5f %5f' % (d,t_FLINT,t_RAZIN,t_FLINT/t_RAZIN)
 else:
  print '%3d: %5f %5f' % (d,t_FLINT,t_RAZIN)

def test_for_dim( dim ):
 global t_FLINT, t_RAZIN
 t_FLINT=t_RAZIN=0
 for i in range(10):
  m=sage.all.random_matrix( sage.all.ZZ, dim, x=-2**1000, y=2*1000 )
  # in 50% cases run flint.det before flint.det_20140704
  test_m(m,i&1)
 print_statistic(dim)
  
def test_m(m,i_m_first):
 global t_FLINT, t_RAZIN
 a=flint.fmpz_mat( m )
 if i_m_first:
  t0=time.time()
  dG = flint.det(a)
  t1=time.time()
  dB = flint.det_20140704( a )
  t2=time.time()
  t_FLINT += t1-t0
  t_RAZIN += t2-t1
 else:
  t0=time.time()
  dB = flint.det_20140704( a )
  t1=time.time()
  dG = flint.det(a)
  t2=time.time()
  t_FLINT += t2-t1
  t_RAZIN += t1-t0
 if dG==dB:
  return
 print 'test failed, det good/bad=%s/%s' % (dG,dB)
 print m.__repr__()
 assert 0

for d in range(155,256):
 test_for_dim( d )

print '\ntest passed'
