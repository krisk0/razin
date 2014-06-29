#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests subroutine nmod_mat_det_mod_pk()
'''

import sage.all
import flint_sage as flint,sys

randint=sage.all.randint

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

def test_p( dim, p ):
 p=flint.prev_prime( p )
 k_max=max_degree( p )
 if k_max>1:
  test_pk(dim, p, 1 )
  if k_max>2:
   test_pk(dim, p, 2 )
   if k_max >= 4:
    test_pk(dim, p, randint(3,k_max-1) )
 test_pk(dim, p, k_max )

def test_pk( dim, p, k ):
 p_deg_k = p**k
 for i in range(100):
  m=sage.all.random_matrix( sage.all.ZZ, dim, x=0, y=p_deg_k )
  dG = m.det() % p_deg_k
  if dim>4 and 0:
   print 'p=%s k=%s p**k=%s' % (p,k,p_deg_k)
   print 'm=\n',m
  dB = flint.det_mod_pk_3arg( flint.fmpz_mat(m), p, k )
  if dG==dB:
   continue
  print 'p=%s k=%s p**k=%s' % (p,k,p_deg_k)
  print 'test failed, det good/bad=%s/%s' % (dG,dB)
  print m
  assert 0

for d in range(1,15):
 sys.stdout.write('%X ' % d)
 test_for_dim( d )

print '\ntest passed'
