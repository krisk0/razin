#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

import sage.all
import flint_sage
Integer=sage.all.Integer

def test(a,b):
 c=a%b
 d=flint_sage.AmodB(Integer(a),Integer(b))
 print '%X %% %X = %X = %X' % (a,b,c,d)
 assert c==d

test(3,2)
test(2**64-3,2**64-5)

print 'test passed'
