# This program is part of RAZIN
# Licence: GNU General Public License (GPL)

cdef extern from 'flint/fmpq_mat.h':
 void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols)
 void fmpq_mat_clear(fmpq_mat_t mat)

def raw_entry( n, d ):
 if d==1:
  return str(n)
 if d<0 and n<0:
  d,n=-d,-n
 return '%s/%s' % (d,n)

cdef class fmpq_mat:

 cdef fmpq_mat_t matQQ

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
 
 def __init__(self, m):
  '''
   m is fmpq_mat
   or tuple
    (q,t) where q is Sage Integer and i is fmpz_mat
    (rows,cols,li) where (rows,cols) are dimensions and li is list or array of
     Sage Rational
     
   contents of raw_str() shows that 2-argument form calls gcd(,q) for all entries
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
    #print '__init__: q,m=%s,\n%s' % (q_sage,m)
    fmpz_init( q )
    fmpz_set_mpz( q, (<Integer>q_sage).value )
    #print '__init__: q set'
    fmpq_mat_init(self.matQQ, (<fmpz_mat>t).matr[0].r, (<fmpz_mat>t).matr[0].c )
    #print '__init__: mat_init done'
    size = (<fmpz_mat>t).matr[0].r * (<fmpz_mat>t).matr[0].c
    #print '__init__: size=%s' % size
    for i in range(size):
     fmpq_set_fmpz_frac( self.matQQ[0].entries+i, (<fmpz_mat>t).matr[0].\
      entries+i, q )
    #print '__init__: array initialized'
    fmpz_clear( q )
   else:
    rows,cols,li=m
    fmpq_mat_init(self.matQQ, <long>rows, <long>cols )
    size = <long>rows * <long>cols
    for i in range(size):
     fmpq_set_mpq( self.matQQ[0].entries+i, (<Rational>li[i]).value )
 
 def __dealloc__(self):
  fmpq_mat_clear(self.matQQ)
