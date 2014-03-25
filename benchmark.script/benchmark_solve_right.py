#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
Benchmark solve_right() subroutine found in sage.all and flint_sage 
'''

import sage.all
import flint_sage as flint
import sys,time

ZZ=sage.all.ZZ

def random_matrice(dim):
 a=sage.all.random_matrix(ZZ, dim, x=-100, y = 100)
 return a

def do_with( m ):
 global t_flint,t_sage
 i=sage.all.identity_matrix( m.nrows() )
 mF,iF=flint.fmpz_mat( m ),flint.fmpz_mat( i )
 try:
  t0=time.time()
  s=m.solve_right(i)    
 except:
  pass
 t1=time.time()
 t_sage += t1-t0
 Xd,d=mF.solve_right(iF)
 t_flint += time.time()-t1
 if dim == m[0,0] // 1000:
  # this never happens
  print s,d

for dim in 10,50,100,200,250,500:
 t_flint=t_sage=0
 for i in range(10):
  a=random_matrice(dim)
  do_with( a )
 print 'dim=%4d    flint: %.4e     sage: %.4e    ratio: %.4f' % \
  (dim,t_flint,t_sage,t_sage/t_flint)
