#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests C subroutine fmpz_mat_det_hermitian_decomposition()

for many matrix, det_20140704() should give same result as 
 det_hermitian_decomposition()
'''

import sage.all
import flint_sage as flint
import sys,time
from test_nmod_HNF import unimodular,unimodular_triU

Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

randint=sage.all.randint
random_matrix,identity_matrix=sage.all.random_matrix,sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat
write=sys.stdout.write

count_det=flint.det_20140704
count_hd_det=flint.det_hermitian_decomposition

bits = 99

def generate_matrice(n):
 m=random_matrix(ZZ,n,x=-bits,y=bits)
 return m,count_det( fmpz_mat(m) )

def generate_singular_matrice(n):
 m=random_matrix(ZZ,n,x=-bits,y=bits)
 m[-1]=m[0]-m[1]
 return m

def det_tri(y):
 x=y.export_sage()
 r,d=x.nrows(),1
 for i in range(r):
  d *= x[i,i]
 return d

def would_be_multiplier(a):
 b=flint.nmod_mat_set_fmpz_mat_mod_thalf(a)
 nmod_mat_HNF(b)
 return b.export_nonnegative_fmpz_mat_upper()

def divide_away(m, i, det_i):
 try:
  mm=m.export_sage()
 except:
  mm=m
 try:
  ii=i.export_sage()
 except:
  ii=i
 r = mm * ii
 n=m.nrows()
 for k in range(n):
  for j in range(n):
   r_kj=r[k,j]
   if r_kj % det_i:
    print 'r=\n',r
    print 'det_hermitian_subr(): r[%s,%s] not multiple of %s' % (k,j,det_i)
    assert 0
   r[k,j] = r_kj / det_i
 return r

def decompose_hermitian(m):
 '''
 return 
  c,d such that det c == 1 mod 2 and m*det c = det m
  or None,d if 2**d divides det m
 '''
 a=fmpz_mat(m)
 b=would_be_multiplier(a)
 d0=det_tri(b)
 i,di=small_det_inv( b )
 aa=divide_away( m, i, di )
 if d0<2**63:
  return aa,d0
 aa=fmpz_mat(aa)
 b=would_be_multiplier( aa )
 d1=det_tri(b)
 if d1<2**63:
  i,di=small_det_inv( b )
  return divide_away( aa, i, di ), d0 * d1
 return None, d0 * d1

def count_log2(x):
 return int(x).bit_length()-1

def det_hermitian_time(x):
 t0=time.time()
 rez=count_hd_det( fmpz_mat(x) )
 return rez, time.time()-t0

def det_hermitian(x):
 return count_hd_det( fmpz_mat(x) )

def matrix_with_det( n, d ):
 x=identity_matrix(n)
 x[0,0]=d
 return x * unimodular(n)

def test_for_dim(n):
 '''
 do with small upper-triangular matrix
 
 generate 10 non-singular matrix and one singular
 for x in range 2..100, generate matrice with determinant 2**x and y*2**x,
  where y is a random number in range 3..33
 
 for all the matrix, count its determinant two ways, compare results
 '''
 print 'dim=%s triU test start' % dim
 for i in range(10):
  test_matrice( unimodular_triU(n,5) )
  for k in range(n):
   a=unimodular_triU(n,5)
   for j in range(k+1):
    a[j,j]=randint(-7,7)
   test_matrice(a)
  
 print 'dim=%s singular-nonsingular test start' % dim
 global t0,t1
 t0=t1=0
 singular_count=i=0
 while 1:
  a,d0=generate_matrice(n)
  if d0==0 and singular_count:
   continue
  if d0==0:
   singular_count += 1
  else:
   i += 1
  if n<11 and 0:
   print 'a=\n',a
   print 'a det=%s=%s' % (d0,sage.all.factor(d0))
  test_matrice_mind_time(a)
  if i==10:
   break
 if singular_count==0:
  a=generate_singular_matrice(n)
  test_matrice_mind_time(a)
 print 'test1 t0,1=',t0,t1
 
 t0=t1=0
 for x in range(2,301):
  test_this_det( n, 1<<x )
  test_this_det( n, randint(3,34)<<x )
 print 'test2 t0,1=',t0,t1
  
def test_this_det( n, y ):
 a=matrix_with_det( n, y )
 test_matrice_mind_time( a )
 
def test_matrice_mind_time(a):
 global t0,t1
 m0=time.time()
 z=count_det(fmpz_mat(a))
 m1=time.time()
 #print 'source matrice A with det=%d:\n' % z,a
 h,m2=det_hermitian_time(a)
 if h != z:
  print 'det mismatch: z=%X != %X' % (z,h)
  assert 0
 t0 += m1-m0
 t1 += m2

def test_matrice(x):
 # print 'source matrice:\n',x
 g=count_det(fmpz_mat(x))
 b=count_hd_det( fmpz_mat(x) )
 if g!=b:
  print 'det mismatch, good/bad=%s/%s' % (g,b)
  assert 0

if 1:
 print 'testing simple diagonal matrice'
 x=identity_matrix(5)
 test_matrice(x)
 for i in 2,3,6:
  x[0,0]=1+2*i
  test_matrice(x)
 for i in range(2,200):
  x[0,0] = 1<<i
  test_matrice(x)

sage.all.set_random_seed('20140831')

t0=t1=0
for dim in range(3,30):
 test_for_dim(dim)

print '\ntest passed'
