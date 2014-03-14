# This program is part of RAZIN
# Licence: GNU General Public License (GPL)

cdef extern from 'flint/fmpq_mat.h':
 void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols)
 void fmpq_mat_clear(fmpq_mat_t mat)
 int fmpq_mat_equal(const fmpq_mat_t a,const fmpq_mat_t b)
 void fmpq_mat_scalar_div_fmpz(fmpq_mat_t tgt,const fmpq_mat_t sou,
  const fmpz_t f)
 int  fmpq_mat_get_fmpz_mat(fmpz_mat_t tgt,const fmpq_mat_t sou)
 void fmpq_mat_set_fmpz_mat(fmpq_mat_t tgt,const fmpz_mat_t sou)
 void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t tgt,const fmpz_mat_t sou,
  const fmpz_t d)
 int fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
 int fmpq_mat_inv(fmpq_mat_t X, const fmpq_mat_t A)

def raw_entry( n, d ):
 if d==1:
  return str(n)
 if d<0 and n<0:
  d,n=-d,-n
 return '%s/%s' % (n,d)

cdef class fmpq_mat:

 cdef fmpq_mat_t matQQ

 def __init__(self, m):
  '''
   m is fmpq_mat
   or tuple
    (q,t) where q is Sage Integer and i is fmpz_mat
    (rows,cols,li) where (rows,cols) are dimensions and li is list or array of
     Sage Rational or list or array of smth that can be converted to rational
     
   contents of raw_str() shows that 2-argument form sets matrice into canonical 
    form
  '''
  cdef Py_ssize_t size,i
  cdef fmpz_t q
  if isinstance(m,fmpq_mat):
   fmpq_mat_init(self.matQQ, (<fmpq_mat>m).matQQ.r, (<fmpq_mat>m).matQQ.c )
   size = (<fmpq_mat>m).matQQ[0].r * (<fmpq_mat>m).matQQ[0].c
   for i in range(size):
    fmpq_set( self.matQQ[0].entries+i, (<fmpq_mat>m).matQQ[0].entries+i )
  else:
   if len(m)==2:
    q_sage,t=m
    fmpq_mat_init(self.matQQ, (<fmpz_mat>t).matr[0].r, (<fmpz_mat>t).matr[0].c )
    if q_sage == 1:
     fmpq_mat_set_fmpz_mat( self.matQQ, (<fmpz_mat>t).matr )
    else:
     fmpz_init( q )
     fmpz_set_mpz( q, (<Integer>q_sage).value )
     fmpq_mat_set_fmpz_mat_div_fmpz( self.matQQ, (<fmpz_mat>t).matr, q )
     fmpz_clear( q )
   else:
    rows,cols,li=m
    fmpq_mat_init(self.matQQ, <long>rows, <long>cols )
    size = <long>rows * <long>cols
    for i in range(size):
     fmpq_set_mpq( self.matQQ[0].entries+i, (<Rational>li[i]).value )
 
 def __repr__(self):
  ' this function is not for speed, but for sanity check '
  return "fmpq_mat(%i, %i, [%s])" % (self.matQQ[0].r, self.matQQ[0].c,
              self.entries_str() )

 def entries_str(self):
  ' this function is not for speed, but for sanity check '
  cdef Py_ssize_t size,i
  cdef char* s
  size=self.matQQ[0].r * self.matQQ[0].c
  for i in range(size):
   s=fmpq_get_str( NULL, 10, self.matQQ[0].entries+i )
   S=s
   S=S.replace(' ','')
   if i==0:
    r = S
   else:
    r += ','+S
  return r
 
 def raw_str(self):
  '''
  extract all numerators and denominators directly
  return result as a space-separated string
  '''
  cdef Py_ssize_t size,i
  cdef Integer num,den
  num,den=Integer(0),Integer(0)
  cdef mpq_t e
  mpq_init( e )
  size=self.matQQ[0].r * self.matQQ[0].c
  for i in range(size):
   fmpq_get_mpq( e, self.matQQ[0].entries+i )
   mpz_set( num.value, mpq_numref( e ) )
   mpz_set( den.value, mpq_denref( e ) )
   s=raw_entry( num, den )
   if i==0:
    r = s
   else:
    r += ' '+s
  mpq_clear( e )
  return r
 
 def __dealloc__(self):
  fmpq_mat_clear(self.matQQ)

 def __richcmp__(fmpq_mat a, fmpq_mat b, int op):
  cdef bint r
  if op != 2 and op != 3:
   raise TypeError('fmpq_mat.__richcmp__(): <> relation not defined')
  r = fmpq_mat_equal( a.matQQ, b.matQQ )
  if op == 3:
   r = not r
  return r

 def solve_dixon(self,fmpq_mat B):
  '''
   if self is singular, returns None

   else returns X such that self * X = B
  '''
  cdef fmpq_mat r=fmpq_mat.__new__(fmpq_mat)
  fmpq_mat_init( r.matQQ, B.matQQ[0].r, B.matQQ[0].c )
  if fmpq_mat_solve_dixon(r.matQQ, self.matQQ, B.matQQ ):
   return r
   
 def inverse(self):
  '''
  return inverse of self, or None
  '''
  cdef fmpq_mat r=fmpq_mat.__new__(fmpq_mat)
  fmpq_mat_init( r.matQQ, self.matQQ[0].r, self.matQQ[0].c )
  if fmpq_mat_inv(r.matQQ, self.matQQ):
   return r

def scalar_div_fmpq_3arg(fmpq_mat tgt, fmpq_mat sou, Integer fA):
 '''
 sets tgt to sou/f
 
 tgt must be already initilized and have same dimesion as sou
 
 test show that this procedure canonicalizes result
 '''
 cdef fmpz_t f
 fmpz_init( f )
 fmpz_set_mpz( f, fA.value )
 fmpq_mat_scalar_div_fmpz( tgt.matQQ, sou.matQQ, f )
 fmpz_clear( f )

