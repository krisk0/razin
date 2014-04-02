#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests fmpz_get_junior_limb() subroutine
'''

import sage.all
import flint_sage as flint

Integer=sage.all.Integer

def test_num(x):
 n=Integer(x)
 r0=n % 2**64
 r1=flint.take_mod_2_64(n)
 print 'Sage : %X' % r0
 print 'flint: %X' % r1
 assert r0 == r1

test_num(-8523048840856612143)
test_num( 8523048840856612143)
test_num(-8422362453803)
test_num( 8422362453803)

assert flint.check_fmpz_to_t()
print 'Test passed'
