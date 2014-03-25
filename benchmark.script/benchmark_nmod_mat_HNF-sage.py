#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks my C procedure nmod_mat_HNF() (wrapped as 
fmpz_mat_hermite_form()) against Sage Cython 
Matrix_integer_dense._hnf_modn_impl (wrapped as Matrix_integer_dense._hnf_mod()
'''

import sage.all
import sys,time
try:
 import flint_sage as flint
except:
 print 'you forgot to install flint_sage Python wrapper'
 sys.exit(1)

import sage.matrix.matrix_integer_dense_hnf
double_det=sage.matrix.matrix_integer_dense_hnf.double_det

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=20
dim_increases=19
volume=10
t_sage=0
t_mine=[0,0,0]
coin_fell_on_the_edge=[0,0]
shut_up=1

def random_data(dim):
 ' returns matrice and absolute value of its determinant '
 global coin_fell_on_the_edge
 n = m = dim
 while 1:
  A=sage.all.random_matrix( ZZ, dim, x=-100, y=100 )
  B = A.matrix_from_rows(range(m-2)).matrix_from_columns(range(n-1))
  c = A.matrix_from_rows([m-2]).matrix_from_columns (range(n-1))
  d = A.matrix_from_rows([m-1]).matrix_from_columns (range(n-1))
  try:
   (d1,d2) = double_det (B,c,d, proof=True)
  except:
   print 'coin fell on the edge, point 0'
   coin_fell_on_the_edge[0] += 1
   continue
  (g,k,l) = d1._xgcd (d2, minimal=True)
  assert g >= 0
  if not shut_up:
   print "Stein double-det g=%s"%g
  CUTOFF = 2**30
  # aint not interesting to count HNF of unimodular matrice, so g>0
  if g==1 or 2*g > CUTOFF:
   if not shut_up:
    print 'bad det, going back'
   continue
  W = B.stack (k*c + l*d)
  not_g=abs( flint.det(flint.fmpz_mat( W )) )
  if not_g != g:
   print 'coin fell on the edge, real det=%d != %d' (not_g,g)
   coin_fell_on_the_edge[1] += 1
  return W,g

def test_with( W, g ):
 global t_sage,t_mine
 t0=time.time()
 sage = W._hnf_mod(2*g)
 t_sage += time.time()-t0
 for i in range(3):
  t1=time.time()
  r=flint.fmpz_mat_hermite_form( flint.fmpz_mat( W ), g*(i+1) )
  t_mine[i] += time.time()-t1
  if( r.export_sage() != sage ):
   print 'result incorrect, g=%s mult=%s\n' % (g,i+1)
   print W.str()
   print '\n',sage.str()
   print '\n',r.export_sage().str()
   sys.exit(1)

def zero_time():
 global t_sage,t_mine
 t_sage=0
 t_mine=[0,0,0]

def dump_time(d):
 print '%3d   %.2e    %.2e     %.2e     %.2e' % (d,t_sage,t_mine[0],t_mine[1],
  t_mine[2])
 sys.stdout.flush()

def coin_fall_on_the_edge():
 for i in range(2):
  if coin_fell_on_the_edge[i]:
   print 'at stage %s coin fell on the edge %s times' % \
    (i,coin_fell_on_the_edge[i])

sage.all.set_random_seed('20140318')
print ' dim    sage         x1           x2           x3' 
for i in range(dim_increases):
 dim=dim0+i*10
 zero_time()
 for j in range(volume):
  W,g=random_data(dim)
  test_with(W,g)
 dump_time(dim)
coin_fall_on_the_edge()
