#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests strong Miller-Rabin test n_is_big_sureprime()
'''

import sage.all
import flint_sage as flint
import sys

d=flint.primes_in_range( 0xFFFFFFFFFFFFFF4D, 2**64-2 )
for i in sorted( d.keys() ):
 sys.stdout.write( '%X ' % i )
print

print '\ntest passed'
