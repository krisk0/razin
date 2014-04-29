#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests Agrawal-Kayal-Saxena primality test originally implemented by 
 袁轶君 and modified by me
'''

import sage.all
import flint_sage as flint
import sys

AKS=flint.n_is_prime_AKS
Integer=sage.all.Integer

def do_test(x,y):
 print 'testing 0x%X' % y
 assert x==AKS( Integer(y) )

do_test(0,9)
do_test(0,15)
do_test(1,3)
do_test(1,7)
do_test(1,53)
do_test(1,101)
do_test(1,103)

for i in 'C5', 'AD', 'A1', '4D':
 p=(((1<<56)-1) << 8) ^ int(i,16)
 do_test( 1, p )
 do_test( 0, p+2 )

do_test(0,2**64-1)

print '\ntest passed'
