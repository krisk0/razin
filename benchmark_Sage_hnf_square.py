#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program re-implements Sage hnf_square() subroutine for the case when
 input matrice is non-singular, keeping all logic and all features, only
 replacing calls to few low-level subroutines and operations:
   a) S.det() replaced with flint.det( ... )
   b) S.solve_right() replaced with F.solve_dixon()
   c) W._hnf_mod(2*g) replaced with fmpz_mat_hermite_form(W,g)
   d) Some * and / operations on matrice are replaced with fmpq_mat
       multiplivation
   e) Extra test for non-singularity is perfomed in
    solve_system_with_difficult_last_row(), so as not to call solve_dixon() on
    singular matrice
This program also benchmarks the new subroutine against hnf_square(), checks
 that results coincide
'''

import sage.all
import sys,time,numpy,re
try:
 import flint_sage as flint
except:
 print 'you forgot to install flint_sage Python wrapper'
 sys.exit(1)

debug_mode=0

write=sys.stdout.write

fmpq_mat,fmpz_mat=flint.fmpq_mat,flint.fmpz_mat
column_to_fmpq_mat=flint.column_to_fmpq_mat

ZZ=sage.all.ZZ
Integer=sage.all.Integer
import sage.matrix.matrix_integer_dense_hnf
det_given_divisor=sage.matrix.matrix_integer_dense_hnf.det_given_divisor
add_row=sage.matrix.matrix_integer_dense_hnf.add_row
add_column=sage.matrix.matrix_integer_dense_hnf.add_column
pivots_of_hnf_matrix=sage.matrix.matrix_integer_dense_hnf.pivots_of_hnf_matrix
matrix=sage.all.matrix
random_matrix=sage.all.random_matrix

t_sage_min=t_sage_max=0
t_mine_min=t_mine_max=0
bits_choice=[8,32,128,256,512]
#           dim  min_bits max_bits
dim_data=[(  50,     8,     512   ),
          ( 100,     8,     512   ),
          ( 250,     8,     512   ),
          ( 500,     8,     256   ),
          (1000,     8,     128   ),
          (2000,     8,       8   ),
          (4000,     8,       8   )]
five=5
if debug_mode:
 dim_data=[(   4,     8,     512   ),
           (   5,     8,     512   ),
           (   6,     8,     512   )]
 five=55

table_to_print=None

def reimplemented_solve_right( A, b ):
 if debug_mode:
  sage_r=A.solve_right(b)
 Af=flint.fmpq_mat( (Integer(1), flint.fmpz_mat(A)) )
 bf=flint.fmpq_mat( (Integer(1), flint.fmpz_mat(b)) )
 mine_r=Af.solve_dixon(bf)
 if debug_mode:
  e=mine_r.export_column().column()
  assert e == sage_r
  print 'solve_right(): check positive'
 return mine_r

def reimplemented_double_det(A, b, c):
 '''
 nearly identical to double_det(...,proof=True), only uses a faster FLINT 
  method instead of sage solve_right()
 '''
 A = A.transpose()
 b = b.transpose()
 c = c.transpose()
 B = A.augment(b)
 v = reimplemented_solve_right( B, -c )
 # ensure that gcd() to canonicalize v and w only runs once
 db = det_given_divisor(B, v.denominator(), proof=True)
 n = A.nrows()-1
 vn = v.export_entry(n,0)
 if vn == 0:
  raise ValueError('coin fell on the edge: v[n-1] is zero')
 vn=1/vn
 w=flint.fmpq_mat_scalar_mul_rational( v, -vn )
   # w = (-1/vn)*v
 w.entry_mul_Rational( n, vn )  # w[n-1] = w[n-1]/vn
 dc = det_given_divisor(A.augment(c), w.denominator(), proof=True)
 if debug_mode:
  assert db==flint.det ( fmpz_mat( B ) )
  assert dc==flint.det ( fmpz_mat( A.augment(c) ) )
  print 'double_det() check positive'
 return db, dc

def reimplemented_solve_system_with_difficult_last_row(B, a):
 ' return result as fmpq_mat '
 C = sage.all.copy(B)
 while 1:
  while 1:
   '''
   replace last row of C with random small numbers
   make sure the new matrice is non-singular
   '''
   if debug_mode:
    print 'solve_system_with_difficult_last_row() makin matrice'
   C[C.nrows()-1] = random_matrix( ZZ, 1, C.ncols() ).row(0)
   if not quick_nonsigular_test( C ):
    continue
   # solve, export, make matrice row, transpose
   x=reimplemented_solve_right( C, a ).export_column().column()
   break
  D = B.matrix_from_rows(range(C.nrows()-1))
  N = D._rational_kernel_iml()
  # original solve_system_with_difficult_last_row() goes into infinite
  #  recursion loop if N.ncols() != 1
  # if this equality ever happens, failed assert is a lot better than infinite
  #  loop
  assert N.ncols() == 1
  k = N.matrix_from_columns([0])
  w = B[-1]
  a_prime = a[-1]
  lhs = w*k
  rhs = a_prime - w * x
  if lhs[0] == 0:
   if debug_mode:
    print 'solve_system_with_difficult_last_row(): coin fell on the edge'
   continue
  break
 alpha = rhs[0] / lhs[0]
 x=x + alpha*k
 if debug_mode:
  assert B*x == a
  print 'solve_system_with_difficult_last_row(): salvation correct'
 return x

def reimplemented_add_column( B, H_B, a ):
 ' H_B is Sage matrice or fmpz_mat '
 z = reimplemented_solve_system_with_difficult_last_row(B, a)
 if hasattr(H_B,'nrows'):
  H_B = fmpz_mat( H_B )
 H_Bf=fmpq_mat( (Integer(1), H_B) )
 H_Bf.mul( column_to_fmpq_mat(z.column(0)) )
 r=H_Bf.export_fmpz_mat().export_sage()
 return r

def reimplement_small_det_HNF(W, g):
 #instead of W._hnf_mod(2*g)
 return flint.fmpz_mat_hermite_form( fmpz_mat( W ), Integer(g) )
 
def reimplemented_hnf_square( A ):
 '''
 A: square non-singular matrice
 
 nearly identical to sage hnf_square()
 '''
 mn = A.nrows()
 B = A.matrix_from_rows(range(mn-2)).matrix_from_columns(range(mn-1))
 c = A.matrix_from_rows([mn-2]).matrix_from_columns (range(mn-1))
 d = A.matrix_from_rows([mn-1]).matrix_from_columns (range(mn-1))
 b = A.matrix_from_columns([mn-1]).matrix_from_rows(range(mn-2))
 # done slicing
 try:
  d1,d2 = reimplemented_double_det( B, c, d )
 except (ValueError, ZeroDivisionError), msg:
  d1 = flint.det ( fmpz_mat_t( B.stack(c) ) )
  d2 = flint.det ( fmpz_mat_t( B.stack(d) ) )
 g,k,l = d1._xgcd (d2, minimal=True)
 W = B.stack (k*c + l*d)
 if g == 0:
  H = W.echelon_form(algorithm='pari')
 else:
  CUTOFF=2**29
  if g>CUTOFF:
   # should never come here
   f = W.gcd()
   g = g / (f**W.nrows())
   if g<=CUTOFF:
     W0 = (W/f).change_ring(ZZ)
     H = reimplement_small_det_HNF(W0, g)
     H *= f
   else:
    #                        ... or here
    raise NotImplementedError("fallback to PARI!")
  else:
   H = reimplement_small_det_HNF(W, g)
 # if H is fmpz_mat, save a penny and let it be
 x = reimplemented_add_column(W, H, b.stack(matrix(1,1,[k*A[mn-2,mn-1] + 
  l*A[mn-1,mn-1]])))
 if debug_mode:
  xsage=add_column(W, H, b.stack(matrix(1,1,[k*A[mn-2,mn-1] + l*A[mn-1,mn-1]])),
   True)
  assert xsage == x
 # if H is fmpz_mat, convert it to Sage
 if not hasattr(H,'nrows'):
  H = H.export_sage()
 Hprime = H.augment(x)
 pivots = range(mn-1)
 if debug_mode:
  assert pivots == pivots_of_hnf_matrix(Hprime)
 Hprime, pivots = add_row(Hprime, A.matrix_from_rows([mn-2]), pivots,
  include_zero_rows=False)
 Hprime, pivots = add_row(Hprime, A.matrix_from_rows([mn-1]), pivots,
  include_zero_rows=False)
 if debug_mode:
  assert Hprime == Hprime.matrix_from_rows(range(mn))
 return Hprime

def do_benchmark( m ):
 '''
 if result mismatches, abort
 
 if time is new record, update t_****_***
 '''
 global t_sage_min,t_sage_max,t_mine_min,t_mine_max
 t0=time.time()
 sage_r=sage.matrix.matrix_integer_dense_hnf.hnf_square(m,True)
 t1=time.time()
 mine_r=reimplemented_hnf_square(m)
 t2=time.time()
 if mine_r != sage_r:
  print 'result incorrect, m=\n',m
  print 'sage result=\n',sage_r
  print '  my result=\n',mine_r
  sys.exit(1)
 t0=t1-t0
 if t0<t_sage_min:
  t_sage_min=t0
 else:
  if t0>t_sage_max:
   t_sage_max=t0
 t1=t2-t1
 if t1<t_mine_min:
  t_mine_min=t1
 else:
  if t1>t_mine_max:
   t_mine_max=t1

def random_data(dim,bits):
 '''
 returns random square non-singular matrice with entries  
 uniformly distributed in the interval −2**bits .. 2**bits
 '''
 while 1:
  a=sage.all.random_matrix( ZZ, dim, x=-1<<bits, y=1<<bits )
  if quick_nonsigular_test( a ):
   return a
  # with non-zero probability matrice a is non-singular, so we do extra
  #  check
  if flint.det(b):
   return a

def quick_nonsigular_test( m ):
 '''
 returns 1 iff determinant of m modulo p0*p1*p1*p3 if non-zero

 where p0...p3 are largest 64-bit primes
 '''
 primes_tail=[ 'C5', 'AD', 'A1', '4D' ]
 for x in primes_tail:
  if flint.det_modulo_prime( m, Integer( (((1<<56)-1) << 8) ^ int(x,16) ) ):
   return 1
 return 0

def benchmark( dim, bits, tries, experiment_no, col_no ):
 '''
 run 2 algorithms multiple times (not more than tries time), store time

 early-abort if accumulated t_mine_max>60
 '''
 global t_sage_min,t_sage_max,t_mine_min,t_mine_max
 t_sage_max=t_mine_max=-1
 t_sage_min=t_mine_min=1e77
 for i in range( tries ):
  m=random_data( dim, bits )
  do_benchmark( m )
  if t_mine_max > 60:
   print 'n=%s bits=%s time=%.1f  benchmarks done=%s, skipping further tries' \
    % (dim,bits,t_mine_max,i+1)
   break
 if t_sage_max==-1:
  t_sage_max=t_sage_min
 if t_mine_max==-1:
  t_mine_max=t_mine_min
 save_time( dim, bits, experiment_no, col_no )

def create_table_to_print( d ):
 ' Initialize table to hold benchmark result '
 i=d[0]
 rows=1+2*len(d)
 cols=2
 for j in bits_choice:
  if j >= i[1] and j <= i[2]:
   cols += 1
 t=numpy.resize( numpy.array( [], dtype=object ), (rows,cols) )
 r0=t[0]
 r0[0],r0[1]='','n'
 col=2
 for j in bits_choice:
  if j >= i[1] and j <= i[2]:
   r0[col]='%s bits' % j
   col += 1
 t[1][0]='Sage'
 t[1+len(d)][0]='re-impl'
 for i in range(rows):
  for j in range(cols):
   if t[i][j] == 0:
    t[i][j] = ''
 return t

def save_time( n, b, exp_no, col ):
 ' fill 2 entries in table_to_print '
 print 'savin time into col=%s, rows %s,%s, min time %s,%s' % \
  ( col, 1+exp_no, 1+len(dim_data)+exp_no, t_sage_min, t_mine_min )
 save_time_subr( t_sage_min, t_sage_max, 1+exp_no, col )
 save_time_subr( t_mine_min, t_mine_max, 1+len(dim_data)+exp_no, col )

def save_time_subr( mIn, mAx, row, col ):
 global table_to_print
 if mAx<1:
  fmt='%.2f'
 else:
  fmt='%.1f'
 if mIn==mAx:
  e=fmt % mIn
 else:
  e,e1=fmt % mIn,fmt % mAx
  if e != e1:
   e=e+'–'+e1
   #e=e+'-'+e1
 table_to_print[ row, col ] = e

sage.all.set_random_seed('20140320')
table_to_print=create_table_to_print( dim_data )
result_row=0
for i in dim_data:
 table_to_print[ 1+              result_row ][ 1 ] = str( i[0] )
 table_to_print[ 1+len(dim_data)+result_row ][ 1 ] = str( i[0] )
 col = 2
 for j in bits_choice:
  if j >= i[1] and j <= i[2]:
   benchmark( i[0], j, five, result_row, col )
   col += 1
   if not debug_mode:
    print 'n=%s bits=%s max time sage/mine=%s/%s' % (i[0],j,t_sage_max,\
     t_mine_max)
    sys.stdout.flush()
 result_row += 1

def pretty_print_result( t, f0 ):
 c=t.shape[1]
 f=numpy.resize( numpy.array( [], dtype=object ), c )
 for i in range(c):
  f[i]=decide_format( t, i, f0 )
 for i in range( t.shape[0] ):
  for j in range( c ):
   s=f[j][0] % t[i][j]
   w=f[j][1]
   if len(s) == w:
    x=w-len( s.lstrip() )
    if x>2:
     y=x>>1
     s=s[y:]+( ' ' * y )
   write(s)
   if j<c-1:
    write(' ')
  write('\n')

def decide_format( t, col, f ):
 req=0
 for i in range(t.shape[0]):
  x=('%'+f) % t[i][col]
  cur=len(x)
  if cur>req:
   req=cur
 return '%'+str(req)+f,req

if not debug_mode:
 pretty_print_result( table_to_print, 's' )
print '\n\nTest passed'
