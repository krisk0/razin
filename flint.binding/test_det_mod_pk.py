#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests subroutine inv_mod_pk()
'''

import sage.all
import flint_sage as flint,sys

Integer=sage.all.Integer
WRITE=sys.stdout.write

sage.all.set_random_seed('20140609')

def test_for_dim( d ):
 test_p( d, 3 )
 test_p( d, 5 )
 test_p( d, 2**32 )
 test_p( d, 2**64-2 )

def max_degree(n):
 k,m=1,n
 assert n>1
 while 1:
  m *= n
  if m>2**64:
   return k
  k += 1
 return k

def test_p( dim, p ):
 p=flint.prev_prime( p )
 k=max_degree( p )
 if k>2:
  test_pk(dim, p, 2 )
 test_pk(dim, p, k )

def test_pk( dim, p, k ):
 p_deg_k = p**k
 for i in range(100):
  m=sage.all.random_matrix( sage.all.ZZ, dim, x=0, y=p_deg_k-1 )
  dG = m.det() % p_deg_k
  dB = flint.det_mod_pk_3arg( flint.fmpz_mat(m), p, k )
  if dG==dB:
   continue
  print 'test failed, det good/bad=%s/%s' % (dG,dB)
  print m
  assert 0

for d in range(1,5):
 test_for_dim( d )

print '\ntest passed'
