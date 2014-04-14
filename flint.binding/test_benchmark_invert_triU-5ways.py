#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks 5 ways to invert an upper-triangular matrice

2 ways from sympy: a.inv('LU') and a.inv()

2 ways from Sage: s.I and 
 s._solve_right_nonsingular_square( identity_matrix(dim), check_rank=False )
 
1 way from FLINT: C subroutine fmpq_mat_inv()
 
Matrice is generated as small-det HNF calculated by test_nmod_HNF_nonsquare.py
'''

import sage.all
import flint_sage as flint
import sys,numpy,time,sympy

Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

randint=sage.all.randint
random_matrix,identity_matrix=sage.all.random_matrix,sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat
write=sys.stdout.write

import sage.matrix.matrix_integer_dense_hnf
hnf=sage.matrix.matrix_integer_dense_hnf.hnf

from test_nmod_HNF import unimodular,left_trans
from test_nmod_HNF_nonsquare import small_nums,three_det_divisors,det_divisor

max_det_divisor=2**64-1

def sympy_0(a):
 d=a.nrows()
 A=sympy.Matrix( d, d, a.list() )
 return A,A.inv()

def sympy_1(a):
 d=a.nrows()
 A=sympy.Matrix( d, d, a.list() )
 return A,A.inv('LU')

def sage_0( A ):
 return A,A.I

def sage_1( A ):
 return A,A._solve_right_nonsingular_square\
  ( identity_matrix(A.nrows()), check_rank=False )

def flint_0( A ):
 Aq=flint.fmpq_mat( (Integer(1), A) )
 return A,Aq.inverse()

def multiply_test_fmpz_fmpq( a, b, dim ):
 '''
 a: fmpz_mat
 b: fmpq_mat
 '''
 r=flint.fmpz_mat_mul_by_fmpq_mat(a,b)
 if not r:
  print 'fmpz_mat_mul_by_fmpq_mat(): non-integral result'
  return 0
 if a.export_sage() != identity_matrix( dim ):
  print 'fmpz_mat_mul_by_fmpq_mat(): multiplication result not identity'
  return 0
 return 1

def test_serie_2(be,A,det_div):
 dim=A.ncols()
 for d in det_div:
  assert 1<d and d <= max_det_divisor
  Ff=flint.fmpz_mat_hermite_form_nonsquare( fmpz_mat(A), d )
  F=Ff.export_sage()
  t0=time.time()
  rY0=sympy_0( F )
  t1=time.time()
  rY1=sympy_1( F )
  t2=time.time()
  rS0=sage_0( F )
  t3=time.time()
  F._clear_cache()
  rS1=sage_1( F )
  t4=time.time()
  rF=flint_0( Ff )
  t5=time.time()
  be[-1] += 1
  be[0] += t1-t0
  be[1] += t2-t1
  be[2] += t3-t2
  be[3] += t4-t3
  be[4] += t5-t4
  Yi,Si=sympy.eye(dim),identity_matrix(dim)
  #print rY0[0],'\n*\n',rY0[1],'\n=\n',rY0[0]*rY0[1],'\n=?=\n',Yi
  assert rY0[0]*rY0[1] == Yi
  assert rY1[0]*rY1[1] == Yi
  assert rS0[0]*rS0[1] == Si
  assert rS1[0]*rS1[1] == Si
  assert multiply_test_fmpz_fmpq( rF[0], rF[1], dim )

def test_serie_1(be,m,c,vol):
 if c<30:
  max_k=c+1
 else:
  max_k=30
 for x in range(vol):
  for k in range(1,max_k):
   assert k <= c
   ' create matrice with abs(det)<2**128 and diagonal with k entries>1 '
   nums=small_nums(k)
   A=identity_matrix(ZZ,c)
   if m>c:
    A=A.stack( matrix( m-c, c ) )
   for i in range(k-1):
    A[i,i]=nums[i]
   A[c-1,c-1]=nums[k-1]
   A *= unimodular(c)
   A=left_trans(A,m)
   d=three_det_divisors( nums )
   test_serie_2(be,A,d)
  # don't start new iteration if too many experiments done 
  if be[-1]>=1000:
   return

def allocate_time_array():
 return numpy.array( [float(0),float(0),float(0),float(0),float(0),0], \
  dtype=object )

def test_serie(c):
 rS=allocate_time_array()
 for i in range(3):
  rC=allocate_time_array()
  test_serie_1(rC,c+i,c,3)
  rS += rC
 return rS,rS[-1]

def show_time( dim, t, v ):
 print '%3X:  %.2e  %.2e  %.2e  %.2e  %.2e  %3X' % (dim,\
  t[0]/v,t[1]/v,t[2]/v,t[3]/v,t[4]/v,v)
 sys.stdout.flush()

if __name__ == "__main__":
 sage.all.set_random_seed('20140413')
 print ' dim      y0        y1        s0        s1        f      vol'
 for i in range(3,26):
  t,v=test_serie(i)
  show_time( i, t, v )

print '\ntest passed'
