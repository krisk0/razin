#!/usr/bin/python2 -B
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests subroutine fmpz_triU_small_det_inverse()/
 fmpz_triU_inverse_smallDet()
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

from test_nmod_HNF import unimodular_triU

def test_with( m ):
 a=fmpz_mat(m)
 good=flint.fmpz_mat_inverse( a )
 baad=flint.fmpz_triU_small_det_inverse( a )
 good_i,baad_i=good[0].export_sage(),baad[0].export_sage()
 if good[1] != baad[1] or good_i != baad_i:
  print 'test failed for a=\n',m
  if good[1] != baad[1]:
   print 'd good/bad=%s/%s' % (good[1],baad[1])
  else:
   print 'd=%s' % good[1]
  print 'good inv=\n',good_i
  if good_i != baad_i:
   print 'baad inv=\n',baad_i
  sys.exit(1)

def test_serie(dim):
 test_with( identity_matrix(dim) )
 for i in range(10):
  a=identity_matrix(dim)
  a[1,1] *= randint(2,99)
  a[dim-1,dim-1] *= randint(2,99)
  test_with( a )
 max_d=int( pow( 2**64-1, 1./dim  ) )
 if 0 and max_d>10:
  max_d=10
 for i in range(10):
  a=unimodular_triU(dim,99)
  for j in range(dim):
   a[j,j]=randint( 2, max_d )
   HNFy_col(a,j)
  test_with(a)

def HNFy_col(m,c):
 for i in range(c):
  m[i,c] %= m[c,c]

if __name__ == "__main__":
 a=matrix(ZZ,4,4,[1,0,0,5, 0,6,0,7, 0,0,1,11, 0,0,0,9] )
 test_with( a )
 a=matrix(ZZ,4,4,[1,0,0,5, 0,6,0,7, 0,0,1,11, 0,0,0,6] )
 test_with( a )
 sage.all.set_random_seed('20140417')
 for i in range(4,16):
  test_serie(i)

print '\ntest passed'
