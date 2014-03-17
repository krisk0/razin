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

def random_matrice():
 return sage.all.random_matrix( ZZ, dim, x=-10, y = 10)

def show(a):
 print 'a=\n',a,'\nhnf(a)=\n',a.hermite_form()

 rank,den,h=fmpz_mat(a).rref()
 print 'h0=',h
 h=h.export_sage()
 print rank,den,'\n',h          # nothing resembling HNF

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

sys.exit(0)
print '\ntest passed'
