#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks my C procedure nmod_mat_HNF() (wrapped as 
 fmpz_mat_hermite_form()) against NTL HNF method
'''

import sage.all
import sys,time
try:
 import flint_sage as flint
except:
 print 'you forgot to install flint_sage Python wrapper'
 sys.exit(1)

from test_nmod_HNF import unimodular,left_trans

ZZ=sage.all.ZZ
write=sys.stdout.write
Integer=sage.all.Integer
dim0=20
dim_increases=19
volume=10
t__ntl=0
t_mine=0
matrix,ntl=sage.all.matrix,sage.all.ntl
randint=sage.all.randint
random_matrix=sage.all.random_matrix

def test(a,loud):
 global t__ntl,t_mine
 b=abs( a.determinant() )
 if b<2 or b>=2**64:
  print 'bad determinant, internal error'
  sys.exit(1)
 t=ntl_HNF_mirror(a)
 if loud:
  print 'a=\n',a
 # benchmark
 t0=time.time()
 nmod_r=flint.fmpz_mat_hermite_form( flint.fmpz_mat( a ), b )
 t1=time.time()
 t=t.HNF(b)
 t2=time.time()
 # convert result
 nmod_r=nmod_r.export_sage()
 if loud:
  print 'hnf(a)=\n',nmod_r
 t = ntl_HNF_unmirror( t )
 if nmod_r != t:
  print 'bad result: ntl\n',t
  sys.exit(1)
 # update time
 t__ntl += t2-t1
 t_mine += t1-t0

def ntl_HNF_mirror( a ):
 n=a.ncols()
 b=matrix(n,n)
 for i in range(n):
  for j in range(n):
   b[j,i]=a[j,n-1-i]
 return ntl.mat_ZZ(n,n,b.list())

def ntl_HNF_unmirror( u ):
 n=u.ncols()
 s=matrix(n,n)
 for i in range(n):
  for j in range(n):
   s[j,i]=u[n-1-j,n-1-i]
 return s

def small_nums(x):
 r=[]
 border=1<<(64//x)
 for i in range(x):
  r.append( randint(2,border) )
 return r

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
   A=sage.all.identity_matrix(ZZ,dim)
   for i in range(k-1):
    A[i,i]=nums[i]
   A[dim-1,dim-1]=nums[k-1]
   A *= unimodular(dim)
   A *= unimodular(dim)
   A=left_trans(A,dim)
   A=left_trans(A,dim)
   test(A,loud)

def dump_time(dim):
 print 'time %s %.2f' % (dim, t__ntl/t_mine )
 sys.stdout.flush()

def test_serie(dim):
 loud=(dim<7)
 if loud:
  vol=7
 else:
  vol=17
 print 'dim=',dim
 test_serie_1(dim,vol,loud) 
 dump_time(dim)

a4=matrix(4,[1,14,8,4,4,7,8,12,2,13,12,0,11,0,15,13])
t4=ntl_HNF_mirror(a4)

if 0:
 print 'sage=\n',a4
 print ' ntl=\n',t4
 print a4.determinant(), t4.determinant()
 print 'sage HNF\n',a4.hermite_form(),' ntl HNF\n',t4.HNF()
 print '\n',ntl_HNF_unmirror( t4.HNF() )
assert a4.hermite_form() == ntl_HNF_unmirror( t4.HNF() )

test(a4,1)
sage.all.set_random_seed('20140319')
for i in range(3,13):
 test_serie(i)
print 'test passed'
