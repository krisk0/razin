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
 return sage.all.random_matrix(ZZ, dim, x=-100, y = 100)

for i in range(10):
 a=random_matrice()
 # both forms of constructor work
 if i&1:
  b=flint.fmpz_mat( a )
 else:
  b=flint.fmpz_mat( (dim, dim, straighten(list(a))) )
 assert a==b.export_sage()
 if i==0:
  print b.export_sage()

print '\ntest passed'
