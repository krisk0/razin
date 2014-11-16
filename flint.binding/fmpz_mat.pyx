# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

cdef extern from 'flint/fmpz_mat.h':
 void fmpz_mat_det(fmpz_t det, const fmpz_mat_t m)
 void fmpz_mat_det_modular(fmpz_t det, const fmpz_mat_t m,int proved)
 void fmpz_mat_print(fmpz_mat_t m)
 void fmpz_mat_clear(fmpz_mat_t m)
 void fmpz_mat_init(fmpz_mat_t m, slong rows, slong cols)
 int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const 
  fmpz_mat_t B)
 int fmpz_mat_equal(const fmpz_mat_t a,const fmpz_mat_t b)
 slong fmpz_mat_rref(fmpz_mat_t res, fmpz_t den, const fmpz_mat_t A)
 slong fmpz_mat_fflu(fmpz_mat_t B, fmpz_t den, slong* perm, const fmpz_mat_t A,
  int rank_check)
 void fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A)
 void fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A)
 int fmpz_mat_inv(fmpz_mat_t tgt, fmpz_t den, const fmpz_mat_t sou)
 void fmpz_mat_init_set(fmpz_mat_t tgt, const fmpz_mat_t src)

cdef extern from 'C/fmpz_mat/window_unsh.c':
 void fmpz_triU_inverse_smallDet(fmpz_mat_t T, fmpz_t d, const fmpz_mat_t S)

cdef extern from 'C/fmpz_mat/det_mod_pk.c':
 mp_limb_t fmpz_mat_det_mod_pk_3arg(fmpz_mat_t M,mp_limb_t p,mp_limb_t k)

cdef extern from 'C/fmpz_mat/det_modular_given_divisor_4block.c':
 void fmpz_mat_det_modular_given_divisor_4block(fmpz_t det,const fmpz_mat_t A,
  fmpz_t divisor)

cdef extern from 'C/fmpz_mat/det_modular_given_divisor_4arg.c':
 void fmpz_mat_det_modular_given_divisor_4arg(mpz_t r, mp_limb_t hb,
  mp_limb_t det_mod_T, fmpz_mat_t A)

cdef extern from 'C/fmpz_mat/det.c':
 void fmpz_mat_det_4block(fmpz_t r, const fmpz_mat_t A)

cdef extern from 'C/fmpz_mat/is_singular.c':
 int fmpz_mat_is_singular(const fmpz_mat_t A)

cdef extern from 'C/fmpz_mat/det_suspected_zero.c':
 void fmpz_mat_det_suspected_zero(mpz_t r,const fmpz_mat_t A,const mpz_t W)

cdef extern from 'C/fmpz_mat/det_modular_given_divisor_8arg.c':
 void fmpz_mat_det_modular_given_divisor_8arg(mpz_t det, nmod_mat_t Amod,
  mpfr_t hb, mpfr_prec_t pr, p_k_pk_t& pp, n_primes_rev_t it,
  mp_limb_t xmod, const fmpz_mat_t A)

cdef extern from 'C/fmpz_mat/scalar_divexact_ui_2arg.c':
 void fmpz_mat_scalar_divexact_ui_2arg(fmpz_mat_t R,mp_limb_t d)

cdef extern from 'C/fmpz_mat/det_hermitian_decomposition.c':
 void fmpz_mat_det_hermitian_decomposition(fmpz_t r, const fmpz_mat_t a)

cdef extern from 'C/fmpz_mat/det_odd.c':
 void fmpz_mat_det_odd(fmpz_t r,const fmpz_mat_t a)

cdef extern from 'C/fmpz_mat/hermitian_decomposition_2.c':
 int fmpz_mat_hermitian_decomposition_2(fmpz_mat_t b,fmpz_t r, const fmpz_mat_t m)

cdef extern from 'C/fmpz_mat/Hadamard_Cramer.c':
 slong hadamard_2arg(mpfr_t b,const fmpz_mat_t m)
 mp_limb_t cramer_rule(const mpfr_t den_bound, 
  mpz_square_mat_t A, mpfr_prec_t pr, slong k)

cdef class fmpz_mat:

 cdef fmpz_mat_t matr
 
 def __init__(self, m):
  '''
   m is Matrix_integer_dense or tuple (rows,cols,li) where 
    li is list or array of sage Integers
  '''
  cdef Py_ssize_t size,i
  if isinstance(m,tuple):
   fmpz_mat_init(self.matr,m[0],m[1])
   size=m[0]*m[1]
   m=m[2]
   for i from 0 <= i < size:
    fmpz_set_mpz(self.matr[0].entries+i, (<Integer>m[i]).value )
  else:
   fmpz_mat_init(self.matr, (<Matrix_integer_dense>m)._nrows, \
    (<Matrix_integer_dense>m)._ncols)
   size=(<Matrix_integer_dense>m)._nrows * (<Matrix_integer_dense>m)._ncols
   for i from 0 <= i < size:
    fmpz_set_mpz(self.matr[0].entries+i, \
     (<Matrix_integer_dense>m)._entries[i] )
 
 def __repr__(self):
   return "fmpz_mat(%i, %i, [%s])" % (self.matr[0].r, self.matr[0].c,
            (", ".join(map(str, self.entries()))))

 def nrows(self):
  return self.matr[0].r

 def ncols(self):
  return self.matr[0].c

 def export_sage(self):
  ' export self as sage Matrix_integer_dense '
  cdef Matrix_integer_dense r=Matrix( self.matr[0].r, self.matr[0].c )
  cdef Py_ssize_t i,j,k=0
  cdef slong* on_row
  for i in range(self.matr[0].r):
   on_row=self.matr[0].rows[i]
   for j in range(self.matr[0].c):
    fmpz_get_mpz( r._entries[k], on_row+j )
    k += 1
  return r

 def entries(self):
  ' return flat list of entries converted to Python int '
  cdef Py_ssize_t i,j
  cdef Integer t
  cdef slong* on_row
  r,t=[],Integer(0)
  for i in range(self.matr[0].r):
   on_row=self.matr[0].rows[i]
   for j in range(self.matr[0].c):
    fmpz_get_mpz( t.value, on_row+j )
    r.append( int(t) )
  return r
  
 def solve_right_slave(self,fmpz_mat B):
  cdef fmpz_mat Xd
  Xd = fmpz_mat.__new__( fmpz_mat )
  fmpz_mat_init( Xd.matr, B.matr[0].r, B.matr[0].c )
  cdef Integer d_sage
  cdef fmpz_t d
  fmpz_init( d )
  d_sage = Integer(0)
  #sig_on()
  if fmpz_mat_solve( Xd.matr, d, self.matr, B.matr ):
   fmpz_get_mpz( d_sage.value, d )
   fmpz_clear( d )
   #sig_off()
   return Xd, d_sage
  fmpz_clear( d )
  #fmpz_mat_clear( Xd.matr )  un-commenting this line results in unhandled crash
  ' is this leaking memory? '
  #sig_off()
  return None,d_sage

 def solve_right(self,B):
  '''
   A=self 
   B is fmpz_mat or smth valid for fmpz_mat constructor
   Solve A*X=B for X. Return pair (Xd,d) where Xd is integer matrice
   If A is singular, return (None,0)
   
   B must have proper dimensions
  '''
  if isinstance(B,fmpz_mat):
   return self.solve_right_slave(B)
  C=fmpz_mat(B)
  return self.solve_right_slave(C)

 def __richcmp__(fmpz_mat a, fmpz_mat b, int op):
  cdef bint r
  if op != 2 and op != 3:
   raise TypeError('fmpz_mat.__richcmp__(): <> relation not defined')
  r = fmpz_mat_equal( a.matr, b.matr )
  if op == 3:
   r = not r
  return r

 def __dealloc__(self):
  fmpz_mat_clear(self.matr)

 def rref(self):
  '''
  explore what fmpz_mat_rref() does
  
  return triple rank,den,matrice
  '''
  cdef fmpz_mat r=fmpz_mat.__new__( fmpz_mat )
  fmpz_mat_init( r.matr, self.matr[0].r, self.matr[0].c )
  cdef fmpz_t den
  fmpz_init( den )
  cdef slong rank=fmpz_mat_rref(r.matr, den, self.matr )
  cdef Integer d=Integer(0)
  fmpz_get_mpz( d.value, den )
  fmpz_clear( den )
  return int(rank),d,r
  
 def fflu(self):
  '''
  return triple (B, den, perm) as calculated by fmpz_mat_fflu()
  '''
  cdef slong* p=<slong*>malloc(self.matr[0].r * sizeof(slong))
  cdef slong i
  for i in range(self.matr[0].r):
   p[i]=i
  cdef fmpz_mat b=fmpz_mat.__new__( fmpz_mat )
  fmpz_mat_init( b.matr, self.matr[0].r, self.matr[0].c )
  cdef fmpz_t d
  fmpz_init( d )
  fmpz_mat_fflu(b.matr, d, p, self.matr, 0)
  cdef Integer d_sage=Integer(0)
  p_list=[]
  for i in range(self.matr[0].r):
   p_list.append( int(p[i]) )
  fmpz_get_mpz( d_sage.value, d )
  fmpz_clear( d )
  free(p)
  return b,d_sage,p_list
  
 def determinant( self ):
  cdef fmpz_t d
  cdef Integer r=Integer(0)
  fmpz_init( d )
  fmpz_mat_det( d, self.matr )
  fmpz_get_mpz( r.value, d )
  fmpz_clear( d )
  return r

 def det_divisor( self ):
  cdef fmpz_t d
  cdef Integer r=Integer(0)
  fmpz_init( d )
  fmpz_mat_det_divisor( d, self.matr )
  fmpz_get_mpz( r.value, d )
  fmpz_clear( d )
  return r
  
 def hadamard_bound( self ):
  cdef fmpz_t d
  cdef Integer r=Integer(0)
  fmpz_init( d )
  fmpz_mat_det_bound( d, self.matr )
  fmpz_get_mpz( r.value, d )
  fmpz_clear( d )
  return r
  
 def solve_dixon( self, fmpz_mat rght ):
  '''
  run Dixon solver to get A**-1 * rght where A=self --- non-singular square

  return result as fmpq_matrice
  
  this subroutine forms fmpq_matrice, later fmpq_mat_solve_dixon() goes
   back to fmpz_matrice
  
  TODO: use fmpz_mat_solve_dixon() directly
  '''
  cdef fmpq_mat A=fmpq_mat.__new__( fmpq_mat )
  fmpq_mat_init(A.matQQ, self.matr[0].r, self.matr[0].r )
  fmpq_mat_set_fmpz_mat( A.matQQ, self.matr )
  # should fail before this point if is self is not square
  cdef fmpq_mat B=fmpq_mat.__new__( fmpq_mat )
  fmpq_mat_init( B.matQQ, rght.matr[0].r, rght.matr[0].c )
  fmpq_mat_set_fmpz_mat( B.matQQ, rght.matr )
  # line below should return None if A is singular
  return A.solve_dixon( B )

def det(fmpz_mat i):
 cdef fmpz_t d
 cdef Integer r=Integer(0)
 fmpz_init( d )
 fmpz_mat_det( d, i.matr )
 fmpz_get_mpz( r.value, d )
 fmpz_clear( d )
 #can return either r or int(r)
 return r

def det_20140704(fmpz_mat i):
 cdef fmpz_t d
 cdef Integer r=Integer(0)
 fmpz_init( d )
 fmpz_mat_det_4block( d, i.matr )
 fmpz_get_mpz( r.value, d )
 fmpz_clear( d )
 return r

def det_hermitian_decomposition(fmpz_mat i):
 cdef fmpz_t d
 cdef Integer r=Integer(0)
 fmpz_init( d )
 fmpz_mat_det_hermitian_decomposition( d, i.matr )
 fmpz_get_mpz( r.value, d )
 fmpz_clear( d )
 return r

def det_modular(fmpz_mat i):
 cdef fmpz_t d
 cdef Integer r=Integer(0)
 fmpz_init( d )
 fmpz_mat_det_modular( d, i.matr, <int>1 )
 fmpz_get_mpz( r.value, d )
 fmpz_clear( d )
 return r

def count_det_suspected_zero(fmpz_mat m, Integer w):
 '''
 det m is known to be a multiple of w and suspected to be zero
 
 count det m
 '''
 cdef Integer r=Integer(0)
 fmpz_mat_det_suspected_zero(r.value, m.matr, w.value)
 return r

cdef wrap_fmpz_mat(fmpz_mat_t a):
 ' shallow constructor '
 cdef fmpz_mat A=fmpz_mat.__new__( fmpz_mat )
 A.matr.entries=a.entries
 A.matr.r=a.r
 A.matr.c=a.c
 A.matr.rows=a.rows
 return A

cdef fmpz_mat_permute( slong* P, fmpz_mat src ):
 '''
 r=copy of src
 apply P to rows of r
 return r
 '''
 cdef fmpz_mat_t tgt
 cdef slong i,j,c=src.matr[0].c
 fmpz_mat_init( tgt, src.matr[0].r, c )
 cdef slong* on_src_row
 cdef slong* on_tgt_row
 for i in range( tgt.r ):
  on_src_row = src.matr[0].rows[ P[i] ]
  on_tgt_row =         tgt.rows[   i  ]
  for j in range(c):
   fmpz_set( on_tgt_row+j, on_src_row+j )
 return wrap_fmpz_mat( tgt )

def fmpz_mat_inverse( fmpz_mat A ):
 '''
 A: square
 returns None,None if A is singular
 returns pair Ainv,den such that (Ainv / den) * A = identity
 '''
 cdef slong t=A.matr[0].r
 cdef fmpz_mat_t Ainv
 fmpz_mat_init( Ainv, t, t )
 cdef fmpz_t den
 fmpz_init(den)
 cdef Integer r=Integer(0)
 t=fmpz_mat_inv( Ainv, den, A.matr )
 if t:
  fmpz_get_mpz( r.value, den )
  fmpz_clear(den)
  return wrap_fmpz_mat(Ainv),r
 fmpz_clear(den)
 fmpz_mat_clear( Ainv )
 return None,None

def fmpz_triU_small_det_inverse( fmpz_mat A ):
 '''
 A: square upper-triangular with positive diagonal
 d=det A is in range 2 .. 2**64-1
 returns pair Ainv,den such that (Ainv / den) * A = identity
 den divides d
 '''
 cdef slong t=A.matr[0].r
 if t < 4:
  return fmpz_mat_inverse( A )
 cdef fmpz_mat_t Ainv
 fmpz_mat_init( Ainv, t, t )
 cdef fmpz_t den
 fmpz_init(den)
 fmpz_triU_inverse_smallDet( Ainv, den, A.matr )
 cdef Integer r=Integer(0)
 fmpz_get_mpz( r.value, den )
 fmpz_clear(den)
 return wrap_fmpz_mat(Ainv),r

def fmpz_mat_copy( fmpz_mat a ):
 cdef fmpz_mat_t c
 fmpz_mat_init_set( c, a.matr )
 return wrap_fmpz_mat(c)

def det_mod_pk_3arg( fmpz_mat a, p, k ):
 return fmpz_mat_det_mod_pk_3arg( a.matr, <mp_limb_t>p, <mp_limb_t>k )

def fmpz_mat_is_singular_wr( fmpz_mat a ):
 return fmpz_mat_is_singular(a.matr)

# include 'fmpz_unimodular.pyx'
include 'fmpz_mat_HNF.pyx'
