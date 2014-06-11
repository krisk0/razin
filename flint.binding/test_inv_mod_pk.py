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
n_preinvert_limb=flint.n_preinvert_limb_wr
inv_6arg=flint.inverse_mod_pk
inv_3arg=flint.inverse_mod_pk_3arg
randint=sage.all.randint
WRITE=sys.stdout.write

def prev_prime(x):
 return Integer(flint.prev_prime(x))

def select_a(p,p_deg_k):
 while 1:
  a=randint(1,p_deg_k)
  if a % p != 0:
   return a

def real_test_3arg(p,k):
 p_deg_k = p**k
 a=select_a( p, p_deg_k )
 q=inv_3arg( a, p, k )
 assert q<p_deg_k
 if q*a % p_deg_k != 1:
  test_failed(a,p,k,q)

def test_3arg(p,k):
 real_test_3arg(p,k)
 if k>1:
  k=randint(1,k-1)
  real_test_3arg(p,k)

def test_6arg(p, k, p_deg_k, p_deg_k_norm, p_deg_k_inv):
 a=select_a( p, p_deg_k )
 q=inv_6arg( a, p, k, p_deg_k, p_deg_k_norm, p_deg_k_inv )
 if q*a % p_deg_k != 1:
  test_failed(a,p,k,q)

def test_failed(a,p,k,r):
 print 'test failed, a p k=%s %s %s %s' % (a,p,k,r)
 sys.exit(1)

def max_degree(n):
 k,m=1,n
 assert n>1
 while 1:
  m *= n
  if m>2**64:
   return k
  k += 1
 return k

def shift_left(x):
 while x<2**63:
  x <<= 1
 return x

def test_serie( rb ):
 p=prev_prime( randint(3,rb) )
 k=max_degree(p)
 p_deg_k=p**k
 p_deg_k_norm=shift_left(p_deg_k)
 p_deg_k_inv=Integer( n_preinvert_limb(p_deg_k_norm) )
 for i in range(12):
  test_3arg(p,k)
  test_6arg(p,k,p_deg_k,p_deg_k_norm,p_deg_k_inv)

def test_number(p):
 WRITE( 'near prime for %X is ' % p )
 sys.stdout.flush()
 p=prev_prime(p)
 print '%X' % p
 k=max_degree(p)
 print 'max degree for %X is %X' % (p,k)
 test_3arg(p,k)

sage.all.set_random_seed('20140607')

test_number(3)
test_number(5)
test_number(2**64-2)

for i in range(10):
 test_serie( 1<<10 )
 test_serie( 1<<20 )
 test_serie( 1<<30 )
 test_serie( 1<<32 )
 test_serie( 1<<50 )

for i in range(100):
 test_serie( 1<<63 )
 test_serie( 1<<64-1 )

print '\ntest passed'
