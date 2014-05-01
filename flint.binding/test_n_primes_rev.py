#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests strong Miller-Rabin test n_is_big_sureprime()
'''

import sage.all
import flint_sage as flint
import sys,numpy,time

count_miine,count_flint=flint.count_primes_in_range,\
 flint.count_primes_in_range_2010
list_miine,list_flint=flint.primes_in_range,flint.primes_in_range_2010
randint=sage.all.randint

d=list_miine( 0xFFFFFFFFFFFFFF4D, 2**64-2 )
assert len(d)==4
d=list_miine( 0xFFFFFFFFFFFFFF4D, 0 )
assert len(d)==4

def print_keys( xx ):
 for x in sorted( xx.keys() ):
  sys.stdout.write( '%X ' % x )
 print

def test_benchm( tt, min_a ):
 a=randint(min_a,2**64-5000-59-1)
 b=a+randint(3,5000)
 b &= 0xFFFFFFFFFFFFFFFE
 t0=time.time()
 c0=count_miine(a,b)
 t1=time.time()
 c1=count_flint(a,b)
 t2=time.time()
 if c0 != c1:
  print 'result mismatch for a,b=%X,%X: %X / %X' % (a,b,c0,c1)
  print_keys( list_miine(a,b) )
  print_keys( list_flint(a,b) )
  if c1<30:
   sys.exit(1)
 tt[0] += t1-t0
 tt[1] += t2-t1

assert 4==count_miine( 0xFFFFFFFFFFFFFF4D, 2**64-2 )
if 0:
 print_keys( list_miine( 0xFFFFFFFFFFFFFF4D, 0xFFFFFFFFFFFFFFAD ) )
 print_keys( list_flint( 0xFFFFFFFFFFFFFF4D, 0xFFFFFFFFFFFFFFAD ) )
assert 3==count_miine( 0xFFFFFFFFFFFFFF4D, 0xFFFFFFFFFFFFFFAD )
assert 3==count_flint( 0xFFFFFFFFFFFFFF4D, 0xFFFFFFFFFFFFFFAD )

sage.all.set_random_seed('20140501')

benchm=numpy.array( [0,0], dtype=float )
for x in range(100):
 test_benchm( benchm, 299210839 )
 sys.stdout.write('.')
print '\ntime: %.2f/%.2f' % (benchm[0],benchm[1])

benchm=numpy.array( [0,0], dtype=float )
for x in range(100):
 test_benchm( benchm, 10**16 )
 sys.stdout.write('.')
print '\ntime: %.2f/%.2f' % (benchm[0],benchm[1])

benchm=numpy.array( [0,0], dtype=float )
for x in range(100):
 test_benchm( benchm, 1<<62 )
 sys.stdout.write('.')
print '\ntime: %.2f/%.2f' % (benchm[0],benchm[1])

print '\ntest passed'
