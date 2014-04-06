#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014

import sage.all
import flint_sage as flint
import sys,numpy

ZZ=sage.all.ZZ
identity_matrix=sage.all.identity_matrix
fmpz_mat=flint.fmpz_mat
matrix=sage.all.matrix
write=sys.stdout.write

w_modulo=2**64

'''
Test tmod_mat object

Create unimodular matrice, lose 1 column, compute its P-L-U decomposition 
 modulo 2**64, check that it is correct
'''

def random_matrice(dim,bits):
 m=sage.all.random_matrix(ZZ, dim, dim, algorithm='unimodular')
 return m.matrix_from_columns( range(dim-1) )

def test_matrice(S):
 m=S.nrows()
 aF_again=fmpz_mat( S )
 PR,LU=flint.tmod_mat_PLU( aF_again )
 if PR == None:
  print_tmod_mat_PLU_error( S, None, 'empty result' )
  sys.exit(1)
 L_sage,U_sage=flint.export_L_sage(LU),flint.export_U_sage(LU)
 tmod_mat_PLU_error=0
 aS_again = flint.tmod_mat_permute_fmpz_mat( PR, fmpz_mat(S) ).export_sage()
 if not ((aS_again - L_sage * U_sage) % w_modulo).is_zero():
  tmod_mat_PLU_error=1
  print_tmod_mat_PLU_error( S, LU, 'P*S != L*U' )
  print ' left part:\n',tmod_pretty_m(aS_again)
  print 'right part:\n',tmod_pretty_m(L_sage * U_sage)
 WL_inv=flint.tmod_mat_solver( PR, LU )
 Wti,Lti=flint.export_Wti( WL_inv ),flint.export_Lti( WL_inv )
 W=U_sage.matrix_from_rows( range(m-1) )
 if (Wti*W.T) % w_modulo != identity_matrix(ZZ,m-1):
  if not tmod_mat_PLU_error:
   tmod_mat_PLU_error=1
   print_tmod_mat_PLU_error( S, LU, 'Wti check failed' )
  print 'Wti=\n',tmod_pretty_m(Wti)
  print 'U=\n',tmod_pretty_m(U_sage)
  print 'W=\n',tmod_pretty_m(W)
 if (Lti*L_sage.T) % w_modulo != identity_matrix(ZZ,m):
  if not tmod_mat_PLU_error:
   tmod_mat_PLU_error=1
   print_tmod_mat_PLU_error( S, LU, 'Lti check failed' )
  print 'solver out=\n',tmod_pretty_m( WL_inv.export_sage() )
  print 'Lti=\n',tmod_pretty_m(Lti)
  print 'L.T=\n',tmod_pretty_m(L_sage.T)
 if tmod_mat_PLU_error:
  sys.exit(1)

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

def tmod_symm_abs(x):
 if x & 0x8000000000000000:
  return w_modulo-x
 return x

print 'startin test'

a32=matrix(ZZ,3,[1,-6,4,-21,7,-48])
PR,LU=flint.tmod_mat_PLU( fmpz_mat(a32) )
LU_inv=flint.tmod_mat_solver(PR, LU)
print tmod_pretty_m( LU.export_sage() ),'\n'
print tmod_pretty_m( LU_inv.export_sage() ),'\n'
print tmod_pretty_m( flint.export_Wti(LU_inv) )
print tmod_pretty_m( flint.export_Lti(LU_inv) )

sage.all.set_random_seed('20140402')
for dim in (20,):
 for i in range(20):
  aS=random_matrice(dim,12)
  if i==1:
   # warp matrice so at least one element is big enough
   if aS[1,0] >= 0:
    aS[1,0] += 3 * w_modulo
   else:
    aS[1,0] -= 5 * w_modulo
  test_matrice(aS)
  if 0:
   print 'so far so good\n\n'
  else:
   write('.'),sys.stdout.flush()
 write('+'),sys.stdout.flush()

print '\ntest passed'
