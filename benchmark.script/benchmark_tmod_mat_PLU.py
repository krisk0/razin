#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys,time

ZZ=sage.all.ZZ
identity_matrix=sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat
matrix=sage.all.matrix
write=sys.stdout.write

w_modulo=2**64

'''
benchmark call sequence tmod_mat_PLU() followed by tmod_mat_solver()

Also print bigger values of diagonal, in compressed form
'''

def random_Utria(dim,bits):
 U=sage.all.random_matrix(ZZ,dim,x=-bits,y=bits)
 for i in range(dim):
  if sage.all.randint(0,1):
   U[i,i]=1
  else:
   U[i,i]=-1
  for j in range(i):
   U[i,j]=0
 return U

def random_Ltria(dim,bits):
 return random_Utria(dim,bits).transpose()

def random_matrice(dim,bits):
 m=identity_matrix(ZZ,dim)
 for i in range(dim):
  if sage.all.randint(0,1):
   m[i,i]=-1
 bits = 1<<bits
 m *= random_Ltria(dim,bits)
 m *= random_Utria(dim,bits)
 m *= random_Ltria(dim,bits)
 m *= random_Utria(dim,bits)
 return m.matrix_from_columns( range(dim-1) )

def test_matrice(S,dim):
 t0=time.time()
 PR,LU=flint.tmod_mat_PLU( fmpz_mat( S ) )
 if PR == None:
  print_tmod_mat_PLU_error( S, None, 'empty result' )
  sys.exit(1)
 WL_inv=flint.tmod_mat_solver( PR, LU )
 t0=time.time()-t0
 print 'time=%s' % t0
 d=flint.agnostic_array_export_big( PR, S.nrows(), S.ncols(), 1<<10 )
 pretty_print_big(dim,d)
 return t0

def pretty_print_big(m, dd):
 if dd==None:
  print 'dim=%s     no big values' % m
  return
 bb=dict()
 for k,v in dd.iteritems():
  l=k.bit_length()
  try:
   bb[l] += v
  except:
   bb[l] = v
 print 'dim=%s  dict size=%s   champions:' % (m,len(dd))
 for k,v in bb.iteritems():
  print '    %5X : %5X' % (k,v)
 print
 sys.stdout.flush()

def print_tmod_mat_PLU_error( S, LU, m ):
 print 'test_matrice(): error in tmod_mat_PLU(): '+m+'\n'
 if LU != None:
  print 'LU as came from flint:\n',tmod_pretty_m(LU.export_sage())
 print '\nsource matrice:\n',tmod_pretty_m(S)

def tmod_pretty_m(x):
 y,m,n=sage.all.copy(x),x.nrows(),x.ncols()
 y = y % w_modulo
 for i in range(m):
  for j in range(n):
   c=y[i,j]
   if c > 0x8000000000000000:
    y[i,j] = c-w_modulo
 return y

def tmod_pretty_num(x):
 if x & 0x8000000000000000:
  return x-w_modulo
 return x

def export_sage(x):
 try:
  return x.export_sage()
 except:
  return x

def apply_Q( m, pq ):
 n,c=m.nrows(),m.ncols()
 q=pq[n:]
 for i in range(c):
  j=q[i]
  if j != i:
   m.swap_columns(i,j)
   q[j]=i

n_invmod=sage.all.inverse_mod

def tmod_pretty_a(a,s):
 r=[]
 for i in range(s):
  r.append( a[i] )
 return r

print 'startin benchmark'

a32=matrix(ZZ,3,[1,-6,4,-21,7,-48])
PR,LU=flint.tmod_mat_PLU( fmpz_mat(a32) )
LU_inv=flint.tmod_mat_solver(PR, LU)
print tmod_pretty_m( LU.export_sage() ),'\n'
print tmod_pretty_m( LU_inv.export_sage() ),'\n'
print tmod_pretty_m( flint.export_Wti(LU_inv) )
print tmod_pretty_m( flint.export_Lti(LU_inv) )

sage.all.set_random_seed('20140402')
tries=5
for dim in range(100,3001,50):
#(20,100,250):
 s=0
 for i in range(tries):
  aS=random_matrice(dim,64)
  s += test_matrice(aS,dim)
 print 'for dim=%s avg time=%.2f' % (dim,s/tries)
 sys.stdout.flush()
