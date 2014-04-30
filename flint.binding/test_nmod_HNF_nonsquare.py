#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests nmod_mat_HNF_nonsquare()
'''

import sage.all
import flint_sage as flint
import sys
Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

randint=sage.all.randint
random_matrix,identity_matrix=sage.all.random_matrix,sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat

import sage.matrix.matrix_integer_dense_hnf

from test_nmod_HNF import unimodular,left_trans
hnf=sage.matrix.matrix_integer_dense_hnf.hnf
write=sys.stdout.write

max_det_divisor=2**64-1

def test_dump(a,h,mu):
 print 'hnf of\n',a,'\n modulo',mu,'equals\n',h

def test_serie_2(A,det_div,loud):
 H=hnf(A,include_zero_rows=0)[0]
 for d in det_div:
  assert 1<d and d <= max_det_divisor
  F=flint.fmpz_mat_hermite_form_nonsquare( fmpz_mat(A), d )
  if loud:
   print '\n\n\n'
   test_dump(A, F.export_sage(), d)
   print '\n\n'
  err=0
  if trace(F) != d:
   print 'Determinant of F is bad: %s != %s' % (trace(F),d)
   err=1
  if not divides( F, H ):
   print 'F does not divide H'
   err=1
  if err: 
   if not loud:
    test_dump(A, F.export_sage(), d)
   print 'H=\n',H
   if A.nrows() == A.ncols():
    x=trace(H)
    if x<2**64:
     Q=flint.fmpz_mat_hermite_form( fmpz_mat(A), x ).export_sage()
     print 'Q=\n',Q
   sys.exit(1)

def trace(x):
 try:
  y=x.export_sage()
 except:
  y=x
 r=y[0,0]
 for i in range(1,y.ncols()):
  r *= y[i,i]
 return r

def divides(A, b):
 a = A.export_sage()
 c = b * a.I
 try:
  c=c.change_ring(ZZ)
  return 1
 except:
  return 0

def det_divisor( nn ):
 s=len( nn )
 i=randint( 0, s-1 )
 r=nn[i]
 i0=i
 while 1:
  assert r <= max_det_divisor
  i = (i+1) % s
  #print 'det_divisor(): r=0x%X' % r
  if i==i0 or 0 == randint(0,2):
   return r    # full round passed, or just randomly exit
  r_new = r*nn[i]
  if r_new > max_det_divisor:
   return r
  r=r_new

def three_det_divisors( li ):
 if len(li)==1:
  return [ Integer( li[0] ) ]
 r,i=dict(),0
 while i<50:
  r[ Integer( det_divisor(li) ) ] = 0
  if len(r) >= 3:
   return r.keys()
  i += 1
 return r.keys()

def small_nums(x):
 r=[]
 border=128//x
 if border>64:
  border=64
 if 0 and border>9:
  border=9
 border=(1<<border)-1
 assert border>2
 for i in range(x):
  r.append( randint(2,border) )
 return r

def test_serie_1(m,c,vol,loud):
 if c<30:
  max_k=c+1
 else:
  max_k=30
 for k in range(1,max_k):
  assert k <= c
  for x in range(vol):
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
   if 1:
    write('.')
   else:
    print 'nums=%s divisors=%s' % (nums,d)
   test_serie_2(A,d,loud)

def test_serie(c):
 loud=(c<7)
 if loud:
  vol=7
 else:
  vol=17
 for i in range(3):
  m=c+i
  test_serie_1(m,c,vol,loud) 

if __name__ == "__main__":
 sage.all.set_random_seed('20140413')
 for i in range(3,13):
  test_serie(i)
 print 'test passed'
