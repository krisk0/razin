#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program benchmarks flint .solve_right() and also NTL and FLINT determinant
 calculation
 
Reports GMP or MPIR version 
'''

import sage.all
import flint_sage as flint
import sys,time,os,ctypes,subprocess

ZZ=sage.all.ZZ
Integer=sage.all.Integer
write=sys.stdout.write
fmpz_mat=flint.fmpz_mat
Popen,PIPE=subprocess.Popen,subprocess.PIPE

dim=90
tries=20
bits=100
mpir_ver_length=15

def random_matrix(dim,bits):
 '''
 returns random square non-singular matrice with entries
 uniformly distributed in the interval −2**bits .. 2**bits
 '''
 m=1<<bits
 while 1:
  a=sage.all.random_matrix( ZZ, dim, x=-m, y=m )
  if not flint.fmpz_mat_is_singular_wr( fmpz_mat(a) ):
   return a

det_mismatch=0

def benchmark(dim,bits,tries):
 r0=r1=r2=0
 global det_mismatch
 for i in range(tries):
  m=random_matrix(dim,bits)
  r=sage.all.identity_matrix( dim ).column(1).column()
  t0=time.time()
  fmpz_mat(m).solve_dixon( fmpz_mat(r) )
  t1=time.time()
  d_flint=fmpz_mat(m).determinant()
  t2=time.time()
  d_ntl  =m.determinant(algorithm='ntl',proof=True)
  t3=time.time()
  r0 += t1-t0
  r1 += t2-t1
  r2 += t3-t2
  if d_flint != d_ntl:
   det_mismatch=1
 return r0/tries,r1/tries,r2/tries

def stretch_string(x, s):
 s=len(x)-s
 if s>0:
  x += ' ' * s
 return x

def crunch_sage_ver( x ):
 x=x.replace(',',' ').split(' ')
 for y in x:
  if y=='':
   continue
  if y[0].isdigit():
   return 'Sage '+y
  
def ntl_mpir_ver( sage_file ):
 '''
 torture .so files until they say which MPIR or GMP they use
 '''
 sage_file=sage_file.replace( 'all.pyc', 'rings/integer.so' )
 if os.path.isfile(sage_file):
  ntl_so = which_so( sage_file, '/libntl-' )
 else:
  print 'Failed to find Sage integer.so. Not Linux?'
  return '???'
 if not os.path.isfile( ntl_so ):
  print 'Failed to find ntl .so. Sage without NTL? Impossible'
  return '???'
 mpir_so=which_so( ntl_so, '/libmpir' )
 if mpir_so == '':
  gmp_so=which_so( ntl_so, '/libgmp' )
  return add_string('GMP ',take_substring( gmp_so, '__gmp_version' ))
 else:
  return add_string('MPIR ',take_substring( mpir_so, '__mpir_version' ))

def which_so( master_so, slave_so ):
 #print 'ldd "'+master_so+'"|fgrep '+slave_so
 try:
  row = Popen('ldd "'+master_so+'"|fgrep '+slave_so,shell=True,stdout=PIPE).\
   communicate()[0]
 except:
  print 'Failed to execute ldd on '+master_so
  return ''
 # I don't know what ldd returns when path to library contains spaces. No 
 #  guesses, the code below will fail in that case
 x=row.find(' => ')
 slave=row[x+4:].lstrip().split(' ')[0]
 if os.path.isfile(slave):
  return slave
 #print 'ldd says >'+row+'<, but >'+slave+'< is not a file'
 return ''

def take_substring( so_name, var_name ):
 try:
  L=ctypes.cdll.LoadLibrary(so_name)
  v=ctypes.c_char_p.in_dll(L,var_name)
  return v.value
 except:
  pass

def add_string(a, b):
 if b==None:
  return '???'
 return a+b
  
sage.all.set_random_seed('20140307')
ver0=stretch_string( flint.GMP_or_MPIR_version(), mpir_ver_length)
ver1=stretch_string( crunch_sage_ver( sage.all.version() ), mpir_ver_length )
ver2=stretch_string( ntl_mpir_ver( sage.all.__file__ ), mpir_ver_length )
 
ver=ver0+'|'+ver1+'|'+ver2
h=(' ' * len(ver)) + ' Dixon solve     FLINT det       NTL det'
t=benchmark(dim,bits,tries)

print h
print '%s     %.2f           %.2f           %.2f' % (ver,t[0],t[1],t[2])

if det_mismatch:
 print 'Test failed'
else:
 print 'Test passed'
