# This program is part of RAZIN
# Licence: GNU General Public License (GPL)

cdef extern from 'flint/fmpq_mat.h':
 void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols)
 void fmpq_mat_clear(fmpq_mat_t mat)

cdef class fmpq_mat:

 cdef fmpq_mat_t matQQ

 def __repr__(self):
   return "fmpq_mat(%i, %i, [%s])" % (self.matQQ[0].r, self.matQQ[0].c,
            (  ", ".join( self.entries() )  )   )

 def entries(self):
  cdef Py_ssize_t size,i
  cdef char* s
  r,first='',1
  size=self.matQQ[0].r * self.matQQ[0].c
  for i in range(size):
   s=fmpq_get_str( NULL, 10, self.matQQ[0].entries+i )
  if first:
   r += s
  else:
   r += ' '+s
  return r
 
 def __init__(self, m):
  '''
   m is fmpq_mat
   or tuple
    (q,t) where q is Sage Integer and i is fmpz_mat
    (rows,cols,li) where (rows,cols) are dimensions and li is list or array of
     Sage Rational
  '''
  cdef Py_ssize_t size,i
  cdef fmpz_t q
  if isinstance(m,fmpq_mat):
   fmpq_mat_init(self.matQQ, (<fmpq_mat>m).matQQ.r, (<fmpq_mat>m).matQQ.c )
   size = (<fmpq_mat>m).matQQ[0].r * (<fmpq_mat>m).matQQ[0].c
   for i in range(size):
    fmpq_set( self.matQQ[0].entries+i, (<fmpq_mat>m).matQQ[0].entries+i )
  else:
   if len(m==2):
    q_sage,t=m
    fmpz_init( q )
    fmpz_set_mpz( q, (<Integer>q_sage).value )
    fmpq_mat_init(self.matQQ, (<fmpz_mat>t).matr.r, (<fmpz_mat>t).matr.c )
    size = (<fmpz_mat>t).matr.r * (<fmpz_mat>t).matr.c
    for i in range(size):
     fmpq_set_fmpz_frac( self.matQQ[0].entries+i, (<fmpz_mat>t).matr[0].\
      entries+i, q )
    fmpz_clear( q )
   else:
    rows,cols,li=m
    fmpq_mat_init(self.matQQ, <long>rows, <long>cols )
    size = <long>rows * <long>cols
    for i in range(size):
     fmpq_set_mpq( self.matQQ[0].entries+i, (<Rational>li[i]).value )
 
 def __dealloc__(self):
  fmpq_mat_clear(self.matQQ)
