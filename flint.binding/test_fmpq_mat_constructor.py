#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program tests if 3 forms of constructor produce same result

Also calls raw_str() to print unfiltered entries
'''

import sage.all
import flint_sage as flint
import sys

dim=4
QQ=sage.all.QQ
ZZ=sage.all.ZZ
Integer=sage.all.Integer
Rational=sage.all.Rational
F=flint.fmpq_mat

def straighten( xx ):
 r=[]
 for i in xx:
  for j in i:
   r.append( Rational(j) )
 return r

def random_matrice():
 return sage.all.random_matrix(ZZ, dim, x=-10, y = 10)

def sage_ZZ_to_flint_QQ( d, a ):
 return F( (Integer(d), flint.fmpz_mat( a ) ) )

for i in range(5):
 a=random_matrice()
 a[0,0]=-6
 print 'Sage a=\n',a
 b0=sage_ZZ_to_flint_QQ( 1, a )
 for i in range(dim):
  j=Integer(i)
  assert( b0.extract_numerator( j ) == a[j,0] )
 b1=F( b0 )
 b2=F( (dim,dim,straighten( list(a) ) ) )
 s0,s1,s2=str(b0),str(b1),str(b2)
 assert s0 == s1
 assert s1 == s2
 assert b0 == b1
 assert b1 == b2
 print s0,'\n'
 print b2.raw_str()

 b2=sage_ZZ_to_flint_QQ(3,a)
 flint.scalar_div_fmpq_3arg( b0, b1, Integer(3) )
 print '\nb2=%s' % b2.raw_str() # -2 at the beginning of the line shows that 
                                #  constructor divides away denominator
 print 'b0=%s' % b0.raw_str()
 assert b0 == b2

print 'test passed'
