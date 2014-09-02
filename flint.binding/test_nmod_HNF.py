#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests fmpz_mat_hermite_form() function, which is a wrapper above
  nmod_mat.hermite_form(), which is wrapper above nmod_mat_HNF()
'''

import sage.all
import flint_sage as flint
import sys
Integer=sage.all.Integer
matrix=sage.all.matrix
ZZ=sage.all.ZZ

#import sage.misc.prandom

randint=sage.all.randint
random_matrix,identity_matrix=sage.all.random_matrix,sage.all.identity_matrix

def test_dump(a,h):
 print 'hnf of\n',a,'\nequals\n',h

def test(a,loud):
 b=abs( a.determinant() )
 if b<2 or b>=2**63:
  print 'bad determinant, internal error'
  sys.exit(1)
 sage_r=a.hermite_form()
 for mult in 1,2:
  if loud:
   print '\n\n\n'
   test_dump(a,sage_r)
   print '\n\n'
  nmod_r=flint.fmpz_mat_hermite_form( flint.fmpz_mat( a ), b*mult )
  if nmod_r != flint.fmpz_mat( sage_r ):
   if not loud:
    test_dump(a,sage_r)
   print 'mult=',mult
   nmod_r=nmod_r.export_sage()
   print 'nmod result:\n',nmod_r
   det_sage=mult_diag(sage_r)
   det_nmod=mult_diag(nmod_r)
   print 'b=%X det_sage=%X det_nmod=%X' % (b,det_sage,det_nmod)
   sys.exit(1)

def mult_diag(m):
 r=1
 for i in range(m.nrows()):
  r *= m[i,i]
 return r

def m(x):
 dim=int( len(x)**.5 )
 return matrix( dim, x )

def test_serie_1(dim,vol,loud):
 if dim<20:
  max_k=dim+1
 else:
  max_k=21
 for k in range(1,max_k):
  assert k <= dim
  for x in range(vol):
   ' create matrice with abs(det)<2**63 and diagonal with k entries>1 '
   nums=small_nums(k)
   A=identity_matrix(ZZ,dim)
   for i in range(k-1):
    A[i,i]=nums[i]
   A[dim-1,dim-1]=nums[k-1]
   A *= unimodular(dim)
   A=left_trans(A,dim)
   test(A,loud)

def small_nums(x):
 r=[]
 border=1<<(63//x)
 for i in range(x):
  r.append( randint(2,border) )
 return r

def left_trans(m,dim):
 u=unimodular(dim)
 return u*m

def unimodular_triL(dim, x):
 r=identity_matrix(dim)
 for i in range(1,dim):
  for j in range(i):
   r[i,j] = randint( -x, x )
 return r

def unimodular_triU(dim, x):
 r=identity_matrix(dim)
 for i in range(dim-1):
  for j in range(1+i,dim):
   r[i,j] = randint( -x, x )
 return r

def randomly_permute_rows( m, dim ):
 max_i=dim-1
 for i in range(max_i):
  j=randint( i, max_i )
  if j != i:
   m.swap_rows(i,j)

def unimodular(dim):
 # random_matrix(ZZ, dim, algorithm='unimodular') is too slow
 r=unimodular_triL(dim,99) * unimodular_triU(dim,99)
 randomly_permute_rows(r,dim)
 return r

def test_serie(dim):
 loud=(dim<7)
 if loud:
  vol=7
 else:
  vol=17
 print 'dim=',dim
 test_serie_1(dim,vol,loud) 

if __name__ == "__main__":
 a2=m( [ 5,5,-7,-8] )
 
 u2=m( [ 1,77,0,-1] ) * m( [ 1,0,-87,1] )
 assert abs(u2.determinant()) == 1
 test( a2, 1 )
 test( u2*a2, 1)
 test( a2*u2, 1)
 
 v2=m( [ 1,67,0,-1] )
 b2=a2 * u2 * v2 * m( [ 3,0,30,36] )
 
 test( b2, 1 )
 test( b2 * m( [-1, 17, 0, 1 ] ), 1 )
 
 sage.all.set_random_seed('20140316')

 for i in range(3,13):
  test_serie(i)
 print 'test passed'
