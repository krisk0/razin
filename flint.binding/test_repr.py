#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program solves a system of linear equations in different ways, thus testing 
 wrappers for fmpz_mat_solve, fmpq_mat_solve_dixon

if watch_canonical=1, then print some information if raw_str() returns different
 information (which can only happen when one of the matrix is not in canonical 
 form)
'''

import sage.all
import flint_sage as flint
import sys

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=15
watch_canonical=1

def random_matrice(dim):
 return sage.all.random_matrix( ZZ, dim, x=-100, y=100 )

a=random_matrice(3)
a=flint.fmpq_mat( (Integer(1), flint.fmpz_mat(a)) )
r=str(a)
print 'r =',r
s='aa=r'
exec s
print 'aa=',aa
