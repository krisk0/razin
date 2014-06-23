#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

'''
This program tests subroutine det_mod_pk_fix_SE_corner()
'''

import sage.all
import flint_sage as flint,sys
sys.dont_write_bytecode = True
from test_inv_mod_pk import prev_prime,max_degree

Integer=sage.all.Integer
randint=sage.all.randint
ZZ=sage.all.ZZ
Zn=sage.all.Integers
WRITE=sys.stdout.write

def prepend_row_and_col( a ):
 r=sage.all.matrix(ZZ,4,1,[2,3,4,5]).augment(a)
 r=sage.all.matrix(ZZ,1,5).stack(r)
 r[0,0]=1
 return r

def do_with_matrice( a ):
 a1=prepend_row_and_col(a)
 for p in 3,7,2**32,2**64-2:
  p=prev_prime(p)
  k_max=max_degree(p)
  if k_max>1:
   test_1( a, a1, p, 1 )
   if k_max > 2:
    test_1( a, a1, p, 2 )
    if k_max >= 4:
     test_1( a, a1, p, randint(3,k_max-1) )
  test_1( a, a1, p, k_max )

def die(x):
 print x
 sys.exit(1)

def test_1_failed( a, i, u, negate_det, mess, m ):
 print 'dim5 test failed for a=\n',a % m
 print '\n%s' % mess
 print '\ni=\n',i
 print '\nu=\n',u
 print '\nnegate_det=',negate_det
 sys.exit(1)

def test_1( a, a1, p, k):
 p_deg_k=p**k
 if 0:
  print 'p=%s k=%s p**k=%s' % (p,k,p_deg_k)
  print 'a=\n',a,'\n'
 rez=flint.test_invert_4x4_corner( flint.fmpz_mat(a1), p, k )
 if rez==None:
  a_det=a.determinant() % p
  if a_det==0:
   return
  print 'a=',a
  die('Empty result from test_invert_4x4_corner()')
 negate_det,det,i,u=rez
 i %= p_deg_k
 u %= p_deg_k
 try:
  iInv=i.change_ring(Zn(p_deg_k)).I
 except:
  test_1_failed( a, i, u, negate_det, 'i not invertible', p_deg_k )
 a_p=a.change_ring(Zn(p_deg_k))
 a_det=a.determinant() % p_deg_k
 det %= p_deg_k
 if negate_det&1:
  det=p_deg_k-det
 if det != a_det:
  print "good/bad determinant: %s / %s" % (a_det,det)
  test_1_failed( a, i, u, negate_det, 'bad determinant', p_deg_k )
 if u[1,0]==2 and u[2,0]==3 and u[3,0]==4 and u[4,0]==5:
  if not iInv==a_p:
   print 'p**k=%s' % p_deg_k
   print 'a modulo p**k=\n',a % p_deg_k
   print "i' modulo p**k=\n",iInv % p_deg_k
   test_1_failed( a, i, u, negate_det, 'bad inverse, kind 0', p_deg_k )
 else:
  a1_p=a1.change_ring(Zn(p_deg_k))
  for v in range(4):
   for w in range(4):
    u[1+v,1+w]=iInv[v,w]
  u_p=u.change_ring(Zn(p_deg_k))
  try:
   a_i=a1_p.I
  except:
   test_1_failed( a, i, u, negate_det, 'a1 is not invertible', p_deg_k )
  v=u_p*a_i
  if not is_permutation(v):
   print 'u=\n',u
   test_1_failed( a, i, u, negate_det, 'a is not a row-swapped u', p_deg_k )

def is_permutation(x):
 d=x.determinant()
 if d != 1 and d+1 != 0:
  return 0
 d=x.nrows()
 for i in range(d):
  for j in range(d):
   if x[i,j]>1:
    return 0
 return 1

def grow_element( a ):
 co=10
 while 1:
  i=randint(0,3)
  j=randint(i,3)
  if a[i,j]==0:
   a[i,j]=13
   return
  co -= 1
  if 0==co:
   assert 0

sage.all.set_random_seed('20140620')
a=sage.all.identity_matrix(ZZ,4)
do_with_matrice( a )
for i in range(5):
 grow_element( a )
 do_with_matrice( a )
for i in range(50):
 a=sage.all.random_matrix( ZZ, 4, x=0, y=99 )
 do_with_matrice( a )

print '\ntest passed'
