#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys

dim=9
ZZ=sage.all.ZZ

def straighten( xx ):
 r=[]
 for i in xx:
  for j in i:
   r.append(j)
 return r

def random_matrice():
 a=sage.all.random_matrix(ZZ, dim, x=-100, y = 100)
 return a

sys.stdout.write('+')
for i in range(10):
 a=random_matrice()
 # both forms of constructors work
 if i&1:
  b=flint.fmpz_mat( a )
 else:
  b=flint.fmpz_mat( (dim, dim, straighten(list(a))) )
 if 0:
  print 'b=',repr(b)
 else:
  sys.stdout.write('.')
 assert flint.det(b)==a.determinant(algorithm='ntl',proof=True)

print '\ntest passed'
