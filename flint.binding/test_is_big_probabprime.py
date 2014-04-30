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

Integer=sage.all.Integer

def do_test(x,y):
 print 'testing 0x%X' % y
 assert x==flint.is_big_sureprime( y )

do_test(0,9)
do_test(0,15)
do_test(0,795265023)

for i in 'C5', 'AD', 'A1', '4D':
 p=(((1<<56)-1) << 8) ^ int(i,16)
 do_test( 1, p )
 do_test( 0, p+2 )

do_test(0,2**64-1)

print '\ntest passed'
