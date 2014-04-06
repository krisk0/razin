#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys

dim=3
ZZ=sage.all.ZZ
fmpz_mat=flint.fmpz_mat
matrix=sage.all.matrix

'''
Explore what rref() and fflu() do

Result: 
 rref() outputs garbage
 fflu() correctly counts determinant and prints smth related to source
  matrice, not clear what
  
 fmpz_mat_fflu is documented incorrectly in 4.2.3 pdf: wrong parameter type
'''

def random_matrice():
 return sage.all.random_matrix( ZZ, dim, x=-10, y = 10)

def show(a):
 print 'a=\n',a,'\nhnf(a)=\n',a.hermite_form()

 if 0:
  rank,den,h=fmpz_mat(a).rref()
  h=h.export_sage()
  print rank,den,'rref=\n',h          # result looks like garbage to me
 
 B,den,perm=fmpz_mat(a).fflu()
 B=B.export_sage()
 print 'den=%s perm=%s B=\n' % (den,perm)
 print B

 #rank,den,h=fmpz_mat(a).rref_fraction_free() no such function
 #h=h.export_sage()
 #print rank,den,'\n',h

a=matrix( 2, [5,5,-7,-8] )
show( a )
show( a.transpose() )
for i in range(10):
 a=random_matrice()
 show(a)
 show( a.transpose() )
 print '\n'

show( matrix( 3, [0, 0, 1, 0, 2, 0, 3, 0, 0] ) )

print '\ntest passed'
