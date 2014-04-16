#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks 3 ways to invert an upper-triangular matrice


1 way from Sage: s.I
 
2 ways from FLINT: subroutines fmpq_mat_inv() and fmpz_mat_inv()
 
Matrice is generated as small-det HNF calculated by test_nmod_HNF_nonsquare.py
'''

import sage.all
import flint_sage as flint
import sys,numpy,time

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
from test_nmod_HNF_nonsquare import small_nums,three_det_divisors,det_divisor,\
 trace

max_det_divisor=2**64-1

def sage_0( A ):
 return A,A.I

def flint_0( A ):
 Aq=flint.fmpq_mat( (Integer(1), A) )
 return A,Aq.inverse()

def flint_1( A ):
 r0,r1=flint.fmpz_mat_inverse(A)
 return A,r0,r1

def multiply_test_fmpz_fmpq( a_arg, b, dim ):
 '''
 a: fmpz_mat
 b: fmpq_mat
 '''
 a=flint.fmpz_mat_copy( a_arg )
 r=flint.fmpz_mat_mul_by_fmpq_mat(a,b)
 if not r:
  print 'fmpz_mat_mul_by_fmpq_mat(): non-integral result'
  return 0
 if a.export_sage() != identity_matrix( dim ):
  print 'fmpz_mat_mul_by_fmpq_mat(): multiplication result not identity'
  return 0
 return 1

def multiply_test_fmpz_fmpz_den( A_Ainv_den, I ):
 '''
 A * A_inv must be equal I * den
 '''
 return A_Ainv_den[0].export_sage() * A_Ainv_den[1].export_sage() == \
        I * A_Ainv_den[2]

def test_serie_2(be,A,det_div):
 dim=A.ncols()
 for d in det_div:
  assert 1<d and d <= max_det_divisor
  Ff=flint.fmpz_mat_hermite_form_nonsquare( fmpz_mat(A), d )
  F=Ff.export_sage()
  t2=time.time()
  rS0=sage_0( F )
  t3=time.time()
  rF0=flint_0( Ff )
  t4=time.time()
  rF1=flint_1( Ff )
  t5=time.time()
  #print 'A sou=\n',F
  #print 'A    =\n',rF1[0].export_sage()
  be[-1] += 1
  be[0] += t3-t2
  be[1] += t4-t3
  be[2] += t5-t4
  Si=identity_matrix(dim)
  assert rS0[0]*rS0[1] == Si
  assert multiply_test_fmpz_fmpq( rF0[0], rF0[1], dim )
  # FLINT documentation on mpz_mat_inv() says den must divide determinant
  assert 0 == trace(F) % rF1[2]
  assert multiply_test_fmpz_fmpz_den( rF1, Si )

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
 return numpy.array( [float(0),float(0),float(0),0], dtype=object )

def test_serie(c):
 rS=allocate_time_array()
 for i in range(3):
  rC=allocate_time_array()
  test_serie_1(rC,c+i,c,3)
  rS += rC
 return rS,rS[-1]

def show_time( dim, t, v ):
 print '%3X:  %.2e  %.2e  %.2e    %3X' % (dim,t[0]/v,t[1]/v,t[2]/v,v)
 sys.stdout.flush()

if __name__ == "__main__":
 sage.all.set_random_seed('20140413')
 print ' dim      s0        f0        f1      vol'
 for i in range(3,102):
  t,v=test_serie(i)
  show_time( i, t, v )

print '\ntest passed'
