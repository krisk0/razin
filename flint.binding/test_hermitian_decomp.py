#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests subroutines nmod_mat_HNF() and fmpz_triU_small_det_inverse()
 for modulo=2**63 by counting integer matrice determinant via hermitian 
 decomposition
'''

import sage.all
import flint_sage as flint
import sys,time

Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

randint=sage.all.randint
random_matrix,identity_matrix=sage.all.random_matrix,sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat
write=sys.stdout.write
small_det_inv=flint.fmpz_triU_small_det_inverse
nmod_mat_HNF=flint.nmod_mat_HNF_wr

count_det=flint.det_20140704
bits = 1<<100

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
    print 'det_hermittian_subr(): r[%s,%s] not multiple of %s' % (k,j,det_i)
    assert 0
   r[k,j] = r_kj / det_i
 return r

def det_hermittian_subr(m, i, det_i):
 r=divide_away(m,i,det_i)
 return count_det( fmpz_mat(r) )

def decompose_hermittian(m):
 '''
 return 
  c,d such that det c == 1 mod 2 and m*det c = det m
  or None,d if d divides det m
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

def det_hermitian(x):
 y,z=decompose_hermittian(x)
 if y==None:
  return count_det( fmpz_mat(x) )
 else:
  return count_det( fmpz_mat(y) ) * z

def test_for_dim(n):
 '''
 generate 10 non-singular matrix and one singular
 
 count its determinant two ways, compare results
 '''
 singular_count=i=0
 while 1:
  a,d0=generate_matrice(n)
  if d0==0 and singular_count:
   continue
  if d0==0:
   singular_count += 1
  else:
   i += 1
  assert det_hermitian(a) == d0
  if i==10:
   break
 if singular_count==0:
  a=generate_singular_matrice(n)
  assert det_hermitian(a) == 0

x=identity_matrix(3)
for i in range(2,100):
 x[0,0] = 1<<i
 d=det_hermitian(x)
 assert x[0,0] == d

sage.all.set_random_seed('20140831')

for dim in range(3,30):
 write( '%s ' % dim )
 test_for_dim(dim)

print '\ntest passed'
