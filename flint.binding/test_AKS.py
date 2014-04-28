#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys

AKS=flint.n_is_prime_AKS
Integer=sage.all.Integer

def do_test(x,y):
 assert x==AKS( Integer(y) )

do_test(0,9)
do_test(0,15)
do_test(1,3)
do_test(1,7)
do_test(1,53)
do_test(1,101)
do_test(1,103)
do_test(0,2**64-1)

print '\ntest passed'
