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
 print 'ups, you forgot to install flint_sage Python wrapper'
 sys.exit(1)

import sage.matrix.matrix_integer_dense_hnf
double_det=sage.matrix.matrix_integer_dense_hnf.double_det

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=20
t_sage=0
t_mine=[0,0,0]


def random_data(dim):
 ' returns matrice and absolute value of its determinant '
 n = m = dim
 while 1:
  A=sage.all.random_matrix( ZZ, dim, x=-100, y=100 )
  B = A.matrix_from_rows(range(m-2)).matrix_from_columns(range(n-1))
  c = A.matrix_from_rows([m-2]).matrix_from_columns (range(n-1))
  d = A.matrix_from_rows([m-1]).matrix_from_columns (range(n-1))
  try:
   (d1,d2) = double_det (B,c,d, proof=True)
  except:
   print 'ups, coin fell on the edge, point 0'
   continue
  (g,k,l) = d1._xgcd (d2, minimal=True)
  assert g >= 0
  print "Stein double-det g=%s"%g
  CUTOFF = 2**30
  # aint not interesting to count HNF of unimodular matrice, so g>0
  if g==1 or 2*g > CUTOFF:
   print 'bad det, going back'
   continue
  W = B.stack (k*c + l*d)
  not_g=abs( flint.det(flint.fmpz_mat( W )) )
  if not_g != g:
   print 'coin fell on the edge, real det=%d != %d' (not_g,g)
  return W,g

def test_with( W, g ):
 global t_sage,t_mine
 t0=time.time()
 sage = W._hnf_mod(2*g)
 t_sage += time.time()-t0
 mine,all_equal=[],1
 for i in range(3):
  t1=time.time()
  r=flint.fmpz_mat_hermite_form( flint.fmpz_mat( W ), g*(i+1) )
  t_mine[i] += time.time()-t1
  mine.append( r.export_sage() )
 three_plus=''
 for i in range(3):
  if sage == mine[i]:
   three_plus += '+'
  else:
   three_plus += '-'
 print three_plus
 if three_plus != '+++':
  print 'result incorrect, g=%s' % g
  print W
  sys.exit(1)

def zero_time():
 global t_sage,t_mine
 t_sage=0
 t_mine=[0,0,0]

def dump_time(d):
 print '%3d   %.2e    %.2e     %.2e     %.2e' % (d,t_sage,t_mine[0],t_mine[1],
  t_mine[2])
 sys.stdout.flush()

print ' dim   sage          x1           x2           x3' 
for i in range(19):
 dim=dim0+i*10
 zero_time()
 for j in range(30):
  W,g=random_data(dim)
  test_with(W,g)
 dump_time(dim)
