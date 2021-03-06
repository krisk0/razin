#!/usr/bin/python2
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# Licence: GNU General Public License (GPL)

'''
This program assists in solving riddle 'why time spent in
 solve_system_with_difficult_last_row() is noticeably larger than sum of times of
 kernel calculation and solve_right()'
'''

#todo: output partially filled table after dim=x is finished

import sage.all
import sys,time,numpy
try:
 import flint_sage as flint
except:
 print 'you forgot to install flint_sage Python wrapper'
 sys.exit(1)

debug_mode=0
check_dixon_time=0
benchmark_early_aborts=0

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
randint=sage.all.randint

t_sage_min=t_sage_max=0
t_mine_min=t_mine_max=0
handicap=0
bits_choice=[8,32,64,96,128,196,256,384,512]
#           dim  min_bits max_bits
dim_data=\
         [
          (  50, 8,         512  ),
          ( 100, 8,         512  ),
          ( 250, 8,         512  ),
          ( 500, 8,         256  ),
          (1000, 8,         196  ),
          (2000, 8,          32  ),
          (3000, 8,           8  )
         ]
five=5
if debug_mode:
 dim_data=[(   4,     8,     512   ),
           (   5,     8,     512   ),
           (   6,     8,     512   )]
 five=55

table_to_print=profile_data=None

def profile( e, t ):
 r=time.time()
 global profile_data
 t=r-t
 try:
  s=profile_data[e]
  profile_data[e] = s+t
 except:
  profile_data[e] = t
 return r

def IML_solve_right( a, b ):
 # this code only runs for big and fat matrice
 return a._solve_right_nonsingular_square(b,check_rank=False)

def IML_or_FLINT( a ):
 '''
 returns one iff IML solve_right is expected to be faster than FLINT
  solve_dixon()
 
 Point of balance ought to be found via a polynom or smth else more beautiful
  than un-contiguos function used below
 
 2000 20->19
 1000 19->17
  500 20-35 -> 19-33
 '''
 if check_dixon_time:
  return 0
 n=a.nrows()-1
 if n < 289:     # for dim<290
  return 0       #  only Dixon plays
 # for bigger dim, call IML if entries of a are big enough
 avg_log2=0
 for i in range(10):
  j,k=randint(0,n-1),randint(0,n)
  avg_log2 += int( a[j,k] ).bit_length()
 if check_dixon_time:
  avg_log2 = avg_log2/10.
  global profile_data
  try:
   profile_data[ 'log2_aIJ' ] += avg_log2
  except:
   profile_data[ 'log2_aIJ' ] = avg_log2
 if n <= 401:
  return avg_log2 >= 800    # use Dixon for 80 or less bits
 if n <= 500:
  return avg_log2 >= 400    # use Dixon for 40 or less bits
 if n < 1000:
  # use straight line passing thru points 500,19 and 1000,17
  return avg_log2 >= (5249-n)/25.
 return avg_log2 >= 190     # for 2000 equilibrium is at 20. So fix it here
                            #  for all n>1000

def reimplemented_solve_right( A, b ):
 if debug_mode:
  sage_r=A.solve_right(b)
 mine_r=fmpz_mat(A).solve_dixon( fmpz_mat(b) )
 if debug_mode:
  e=mine_r.export_column().column()
  assert e == sage_r
  print 'solve_right(): check positive'
 return mine_r

def IML_solve_right_or_FLINT_dixon( A, b ):
 ' returns result as sage matrice and time wasted which should be handicaped '
 if check_dixon_time:
  global handicap
  tD = time.time()
  x = IML_solve_right( A, b )
  tE=profile( 'IML solver', tD )-tD
  handicap += tE
  print 'increased handicap by %.2f, new value %.2f' % (tE,handicap)
 t0=time.time()
 if IML_or_FLINT( A ):
  # no double accounting!
  assert not check_dixon_time
  mine_r=IML_solve_right( A, b )
  profile( 'IML solver' , t0 )
 else:
  mine_r=fmpz_mat(A).solve_dixon(fmpz_mat(b)).export_column().column()
  profile( 'dixon-last_row' , t0 )
 if check_dixon_time:
  assert mine_r == x
  return mine_r,tE
 return mine_r,0

def reimplemented_double_det(A, b, c):
 '''
 nearly identical to double_det(...,proof=True), only uses a faster FLINT
  method instead of sage solve_right()
 '''
 t0=time.time()
 A = A.transpose()
 b = b.transpose()
 c = c.transpose()
 B = A.augment(b)
 t1=time.time()
 v = reimplemented_solve_right( B, -c )
 profile( 'solve_right-dd' , t1 )
 # ensure that gcd() to canonicalize v and w only runs once
 t1=time.time()
 db = det_given_divisor(B, v.denominator(), proof=True)
 profile( 'det_given_divisor' , t1 )
 n = A.nrows()-1
 vn = v.export_entry(n,0)
 if vn == 0:
  raise ValueError('coin fell on the edge: v[n-1] is zero')
 vn=1/vn
 w=flint.fmpq_mat_scalar_mul_rational( v, -vn )   # w = (-1/vn)*v
 w.entry_mul_Rational( n, vn )  # w[n-1] = w[n-1]/vn
 t1=time.time()
 dc = det_given_divisor(A.augment(c), w.denominator(), proof=True)
 profile( 'det_given_divisor' , t1 )
 if debug_mode:
  assert db==fmpz_mat( B ).determinant()
  assert dc==fmpz_mat( A.augment(c) ).determinant()
  print 'double_det() check positive'
 profile( 'double_det' , t0 )
 return db, dc

def reimplemented_solve_system_with_difficult_last_row(B, a):
 t0=time.time()
 C = sage.all.copy(B)
 D = B.matrix_from_rows(range(C.nrows()-1))
 tK=profile( 'sswdlr 0', t0 )
 N = D._rational_kernel_iml()
 t1=profile( 'kernel', tK )
 # original solve_system_with_difficult_last_row() goes into infinite
 #  recursion loop if N.ncols() != 1
 # however N.ncols() always equals 1 unless IML counts kernel incorrectly
 # We don't check kernel size
 k = N.matrix_from_columns([0])
 w = B[-1]
 a_prime = a[-1][0]
 lhs = (w*k)[0]
 #assert lhs # .pdf explains that lhs != 0, no need to check
 while 1:
  '''
  replace last row of C with random small numbers
  make sure the new matrice is non-singular
  '''
  if debug_mode:
   print 'solve_system_with_difficult_last_row() makin matrice'
  C[C.nrows()-1] = random_matrix( ZZ, 1, C.ncols() ).row(0)
  if quick_nonsigular_test( C ):
   break
 # solve, export as m*1 matrice
 profile( 'sswdlr 1', t1 )
 x,tE=IML_solve_right_or_FLINT_dixon( C, a )
 tL=time.time()
 rhs = a_prime - (w * x)[0]
 alpha = rhs / lhs
 x=x + alpha*k
 if debug_mode:
  assert B*x == a
  print 'solve_system_with_difficult_last_row(): salvation correct'
 if check_dixon_time:
  t0 += tE
 profile( 'sswdlr 2', tL )
 profile( 'system_with_difficult_last_row' , t0 )
 return x

def reimplemented_add_column( B, H_B, a ):
 ' H_B is Sage matrice or fmpz_mat '
 z = reimplemented_solve_system_with_difficult_last_row(B, a)
 if isinstance(H_B,sage.matrix.matrix_integer_dense.Matrix_integer_dense):
  H_B = fmpz_mat( H_B )
 H_Bf=fmpq_mat( (Integer(1), H_B) )
 H_Bf.mul( column_to_fmpq_mat(z.column(0)) )
 r=H_Bf.export_fmpz_mat().export_sage()
 return r

def reimplement_small_det_HNF(W, g):
 # instead of W._hnf_mod(2*g)
 t0=time.time()
 if g==1:
  r=sage.all.identity_matrix(ZZ, W.nrows())
 else:
  r=flint.fmpz_mat_hermite_form( fmpz_mat( W ), Integer(g) )
 profile( 'small_det_HNF' , t0 )
 return r
 
def reimplemented_hnf_square( A ):
 '''
 A: square non-singular matrice
 
 nearly identical to sage hnf_square()
 '''
 ts=time.time()
 mn = A.nrows()
 B = A.matrix_from_rows(range(mn-2)).matrix_from_columns(range(mn-1))
 c = A.matrix_from_rows([mn-2]).matrix_from_columns (range(mn-1))
 d = A.matrix_from_rows([mn-1]).matrix_from_columns (range(mn-1))
 b = A.matrix_from_columns([mn-1]).matrix_from_rows(range(mn-2))
 # done slicing
 try:
  d1,d2 = reimplemented_double_det( B, c, d )
 except (ValueError, ZeroDivisionError), msg:
  d1 = fmpz_mat_t( B.stack(c) ).determinant()
  d2 = fmpz_mat_t( B.stack(d) ).determinant()
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
 if not isinstance(H,\
         sage.matrix.matrix_integer_dense.Matrix_integer_dense):
  H = H.export_sage()
 Hprime = H.augment(x)
 pivots = range(mn-1)
 if debug_mode:
  assert pivots == pivots_of_hnf_matrix(Hprime)
 t0=time.time()
 Hprime, pivots = add_row(Hprime, A.matrix_from_rows([mn-2]), pivots,
  include_zero_rows=False)
 t1=profile( 'add_row-0', t0 )
 Hprime, pivots = add_row(Hprime, A.matrix_from_rows([mn-1]), pivots,
  include_zero_rows=False)
 profile( 'add_row-1', t1 )
 if debug_mode:
  assert Hprime == Hprime.matrix_from_rows(range(mn))
 profile( 'total', ts )
 return Hprime

def record_bug( o, d ):
 import tempfile
 r=tempfile.NamedTemporaryFile( prefix='/tmp/'+d, delete=0 )
 r.write( '%s' % o.list() )
 close(r)

def do_benchmark( m ):
 '''
 if result mismatches, abort

 return 1+Sage time iff benchmark happened, return 0 if Sage hnf_square() 
  failed
 
 if time is new record, update t_****_***
 '''
 global t_sage_min,t_sage_max,t_mine_min,t_mine_max
 t0=time.time()
 try:
  sage_r=sage.matrix.matrix_integer_dense_hnf.hnf_square(m,True)
 except:
  '''
  Sage was seen failing on line 'x[i,0] = x[i,0]/d' with message
  
  TypeError: no conversion of this rational to integer
  
  Bad result from a subroutine? Maybe H_B is incorrect?
  
  We just get away and create another random matrice
  '''
  print 'Feature in Sage hnf_square(), recording matrice'
  record_bug( m, 'hnf_square-failed-' )
  return 0
 t1=time.time()
 m._clear_cache()
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
  if t_sage_max < 0:
   t_sage_max=t_sage_min
 else:
  if t0>t_sage_max:
   t_sage_max=t0
 t1=t2-t1
 if t1<t_mine_min:
  t_mine_min=t1
  if t_mine_max < 0:
   t_mine_max=t_mine_min
 else:
  if t1>t_mine_max:
   t_mine_max=t1
 return 1+t0

def random_data(dim,bits):
 '''
 returns random square non-singular matrice with entries
 uniformly distributed in the interval −2**bits .. 2**bits
 '''
 m=1<<bits
 while 1:
  a=sage.all.random_matrix( ZZ, dim, x=-m, y=m )
  if not flint.fmpz_mat_is_singular_wr( fmpz_mat(a) ):
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

def find_col( x ):
 for i in range(len(bits_choice)):
  if bits_choice[i] == x:
   return 2+i

def benchmark( dim, bits, tries, experiment_no ):
 '''
 run 2 algorithms multiple times (not more than tries time), store time

 early-abort after 2 or more tries if t_mine_max>60
 '''
 global t_sage_min,t_sage_max,t_mine_min,t_mine_max,profile_data,handicap
 handicap=0
 t_sage_max=t_mine_max=-1
 t_sage_min=t_mine_min=1e77
 avg_sage_time=0
 profile_data=dict()
 col_no=find_col(bits)
 for i in range( tries ):
  while 1:
   m=random_data( dim, bits )
   t_sage=do_benchmark( m )
   if t_sage:
    avg_sage_time += t_sage-1
    break
  if benchmark_early_aborts and t_mine_max > 60 and i:
   tries=i+1
   print 'n=%s bits=%s time=%.1f  benchmarks done=%s, skipping further tries' \
    % (dim,bits,t_mine_max,tries)
   break
 t_mine_min -= handicap/tries  # this is slightly incorrect
 t_mine_max -= handicap/tries  # profile_data['total']-handicap is correct
 profile_output( dim, bits, tries )
 save_time( dim, bits, experiment_no, col_no )
 print 'n=%s bits=%s Sage time=%.2f' % (dim,bits,avg_sage_time)

def sort_in_increasin_order( p ):
 cmp_2nd=lambda x,y: cmp(x[1], y[1])
 l=p.items()
 l.sort( cmp_2nd )
 return l

def detect_profile_format( d ):
 k=0
 for i in d.keys():
  c=len(i)
  if c>k:
   k=c
 return '%'+str(k)+'s %.4f'

def profile_output( dim, bits, tries ):
 global profile_data
 if handicap:
  print 'subtracting from total time %.2f handicap %.2f' % \
   (profile_data['total'],handicap)
  profile_data['total'] -= handicap
 fmt=detect_profile_format( profile_data )
 li=sort_in_increasin_order( profile_data )
 he='********* profile for n=%s bits=%s tries=%s *********' % (dim,bits,tries)
 print '\n\n'+he
 for i in li:
  print fmt % (i[0],i[1])
 print ('*' * len(he))+'\n'
 sys.stdout.flush()

def create_table_to_print( d ):
 ' Initialize table to hold benchmark result '
 i=d[0]
 rows=1+2*len(d)
 cols=2+len(bits_choice)
 t=numpy.resize( numpy.array( [], dtype=object ), (rows,cols) )
 r0=t[0]
 r0[0],r0[1]='','n'
 col=2
 for j in bits_choice:
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
   e=e+'-'+e1
 table_to_print[ row, col ] = e

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
    if x >= 2:
     y=x>>1
     s=s[y:]+( ' ' * y )
     assert len(s)==w
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

if __name__ == "__main__":
 sage.all.set_random_seed('20140327')
 table_to_print=create_table_to_print( dim_data )
 result_row=0
 for i in dim_data:
  table_to_print[ 1+              result_row ][ 1 ] = str( i[0] )
  table_to_print[ 1+len(dim_data)+result_row ][ 1 ] = str( i[0] )
  for j in bits_choice:
   if j >= i[1] and j <= i[2]:
    benchmark( i[0], j, five, result_row )
    if not debug_mode:
     print 'n=%s bits=%s max time sage/mine=%s/%s' % (i[0],j,t_sage_max,\
      t_mine_max)
     sys.stdout.flush()
  result_row += 1
 
 if not debug_mode:
  pretty_print_result( table_to_print, 's' )
 print '\n\nTest passed'
