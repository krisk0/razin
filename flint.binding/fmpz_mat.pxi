# This program is part of RAZIN
# Licence: GPL

cdef extern from 'flint/fmpz_mat.h':
 #long* fmpz_mat_entry(fmpz_mat_t mat, long i, long j)
 void fmpz_mat_det(fmpz_t det, const fmpz_mat_t m)
 void fmpz_mat_print(fmpz_mat_t m)
 void fmpz_mat_clear(fmpz_mat_t m)
 void fmpz_mat_init(fmpz_mat_t m, long rows, long cols)
 int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const 
  fmpz_mat_t B)
 int fmpz_mat_equal(const fmpz_mat_t a,const fmpz_mat_t b)

cdef class fmpz_mat:

 cdef fmpz_mat_t matr
 
 def __init__(self, m): 
  # m is Matrix_integer_dense or tuple (rows,cols,li) where 
  #  li is list or array of sage Integers
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
  #fmpz_mat_print(self.matr)
 
 def __repr__(self):
   return "fmpz_mat(%i, %i, [%s])" % (self.matr[0].r, self.matr[0].c,
            (", ".join(map(str, self.entries()))))

 def entries(self):
  cdef Py_ssize_t size,i
  cdef Integer t
  size=self.matr[0].r * self.matr[0].c
  r,t=[],Integer(0)
  for i in range(size):
   fmpz_get_mpz( t.value, self.matr[0].entries+i )
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
  sig_on()
  if fmpz_mat_solve( Xd.matr, d, self.matr, B.matr ):
   fmpz_get_mpz( d_sage.value, d )
   fmpz_clear( d )
   sig_off()
   return Xd, d_sage
  fmpz_clear( d )
  #fmpz_mat_clear( Xd.matr )  un-commenting this line results in unhandled crash
  sig_off()
  return None,d_sage

 def solve_right(self,B):
  '''
   A=self 
   B is fmpz_mat or smth valid for fmpz_mat constructor
   Solve A*X=B for X. Return pair (Xd,d) where Xd is integer matrice
   If A is singular, return (None,0)
   
   B must have proper dimensions
   
   this
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

def det(fmpz_mat i):
 cdef fmpz_t d
 cdef Integer r
 fmpz_init( d )
 fmpz_mat_det( d, i.matr )
 r=Integer(0)
 fmpz_get_mpz( r.value, d )
 fmpz_clear( d )
 #can return either r or int(r)
 return r
