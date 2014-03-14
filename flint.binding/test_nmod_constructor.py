#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

import sage.all
import flint_sage as flint
import random
Integer=sage.all.Integer
ZZ=sage.all.ZZ

def test(a,b,loud):
 assert b>1
 assert b<2**64
 c=a % b
 n=flint.nmod_mat(flint.fmpz_mat(a),Integer(b))
 d=n.export_fmpz_mat()
 assert d==flint.fmpz_mat(c)
 if loud:
  print c

test( sage.all.matrix( 2, [ 5,5,-7,-8] ), 3, 1 )
x=y=2**64-1
z=2**64

def test_serie(dim,vol,loud):
 for i in range(vol):
  test ( sage.all.random_matrix(ZZ, dim, x=x, y=y), random.randint(2,z), loud )

test_serie(3,3,1)
test_serie(33,333,0)

print 'test passed'
