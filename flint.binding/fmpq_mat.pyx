# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

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
 void fmpq_mat_scalar_div_fmpz(fmpq_mat_t tgt, const fmpq_mat_t sou, 
  const fmpz_t x)
 void fmpq_mat_set(fmpq_mat_t tgt, const fmpq_mat_t sou)
 void fmpq_mat_mul(fmpq_mat_t tgt, const fmpq_mat_t A, const fmpq_mat_t B)
 void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t n, fmpz_t d, const fmpq_mat_t s)
 void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t tgt, const fmpz_mat_t A, 
  const fmpq_mat_t B)
 int fmpq_mat_get_fmpz_mat(fmpz_mat_t tgt, const fmpq_mat_t sou)

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
    (rows,cols,li) where (rows,cols) are dimensions and li is a list or array of
     Sage Rational 
     
   WARNING
    this constructor seg-faults on m not matching the above specification.
    This is a feature, not a bug, 
                                  I WON'T CHANGE THAT
  '''
  cdef Py_ssize_t size,i
  cdef fmpz_t q
  if isinstance(m,fmpq_mat):
   fmpq_mat_init( self.matQQ, (<fmpq_mat>m).matQQ.r, (<fmpq_mat>m).matQQ.c )
   fmpq_mat_set( self.matQQ, (<fmpq_mat>m).matQQ )
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
    try:
     if size>len(li):
      size=len(li)
    except:
     pass
    for i in range(size):
     fmpq_set_mpq( self.matQQ[0].entries+i, (<Rational>li[i]).value )

 def mul( self, fmpq_mat c ):
  '''
  multiplies self by c on the right: a = a * c
  '''
  cdef fmpq_mat_t s
  fmpq_mat_init( s, self.matQQ[0].r, c.matQQ[0].c )
  fmpq_mat_mul( s, self.matQQ, c.matQQ )
  fmpq_mat_clear( self.matQQ )
  self.matQQ[0].entries=s.entries
  self.matQQ[0].r=s.r
  self.matQQ[0].c=s.c
  self.matQQ[0].rows=s.rows

 def export_fmpz_mat( self ):
  '''
  if a=self is integral, export it as fmpz_mat

  else return None
  '''
  cdef fmpz_mat i=fmpz_mat.__new__(fmpz_mat)  
  fmpz_mat_init( i.matr, self.matQQ[0].r, self.matQQ[0].c )
  cdef long ok=fmpq_mat_get_fmpz_mat( i.matr, self.matQQ )
  if ok:
   return i
  print 'non-integral matrice:',self

 def integrify( self ):
  '''
  export self=A as pair (M,d) where M=A*d is fmpz_matrix and d is minimal 
   integer such that M is integer
  '''
  cdef fmpz_mat M=fmpz_mat.__new__(fmpz_mat)  
  fmpz_mat_init( M.matr, self.matQQ[0].r, self.matQQ[0].c )
  cdef fmpz_t den
  fmpz_init( den )
  fmpq_mat_get_fmpz_mat_matwise( M.matr, den, self.matQQ )
  cdef Integer d=Integer(0)
  fmpz_get_mpz( d.value, den )
  return M,d

 def export_column( self ):
  ' export left-most column of self as Sage Vector_rational_dense '
  # TODO: the line below calls Python. How to stay inside Cython?
  cdef Vector_rational_dense r=vector( QQ, self.matQQ[0].r )
  cdef long i
  for i in range(self.matQQ[0].r):
   fmpq_get_mpq( r._entries[i], self.matQQ[0].rows[i] )
  return r

 def export_entry( self, i, j ):
  ' return a[i,j] as Sage Rational '
  cdef Rational e=Rational(0)
  #fmpq_get_mpq( e.value, self.matQQ[0].rows[mpz_get_ui(i.value)]+
  # mpz_get_ui(j.value))
  cdef long ii=<long>i
  cdef long jj=<long>j
  fmpq_get_mpq( e.value, self.matQQ[0].rows[ii]+jj )
  return e

 def extract_numerator( self, Integer i ):
  ' returns numerator of a[i,0] as Sage Integer '
  cdef Integer r=Integer(0)
  #fmpz_get_mpz( r.value, self.matQQ[0].rows[i][0].num )
  cdef fmpq* on_row=self.matQQ[0].rows[i]
  fmpz_get_mpz( r.value, <long*>on_row ) # dirty trick to get .num field
                                         # should be a better way to do this
  return r
 
 def entry_div_Integer( self, i, Integer d ):
  ' changes self: a[i,0] /= d '
  cdef long ii=<long>i
  cdef fmpq* on_row=self.matQQ[0].rows[ii]
  cdef fmpz_t u
  fmpz_init( u )
  fmpz_set_mpz( u, d.value )
  fmpq_div_fmpz( on_row, on_row, u ) # official method called, matrice stays 
                                     #  canonical
  fmpz_clear( u )
 
 def entry_mul_Rational( self, i, Rational m):
  ' changes self: a[i,0] *= m '
  cdef long ii=<long>i
  cdef fmpq* on_row=self.matQQ[0].rows[ii]
  cdef fmpq_t q
  fmpq_init( q )
  fmpq_set_mpq( q, m.value )
  fmpq_mul( on_row, on_row, q )
  fmpq_clear( q )
 
 def denominator( self ):
  ' returns common denominator, assumes self in canonical form '
  cdef fmpz_t lcm
  fmpz_init(lcm)
  fmpz_one(lcm)
  cdef long i,j
  cdef fmpq* on_row
  cdef long c=self.matQQ[0].c
  for i in range(self.matQQ[0].r):
   on_row=self.matQQ[0].rows[i]
   for j in range(c):
    fmpz_lcm(lcm, lcm, (<long*>(on_row+j))+1 ) # how to do this cleanly?
  cdef Integer r=Integer(0)
  fmpz_get_mpz( r.value, lcm )
  fmpz_clear(lcm)
  return r
 
 def __repr__(self):
  ' this function is not for speed, but for sanity check '
  return "fmpq_mat(%i, %i, [%s])" % (self.matQQ[0].r, self.matQQ[0].c,
              self.entries_str() )

 def entries_str(self):
  ' this function is not for speed, but for sanity check '
  cdef Py_ssize_t i,j
  cdef char* s
  cdef fmpq* on_row
  for i in range(self.matQQ[0].r):
   on_row=self.matQQ[0].rows[i]
   for j in range(self.matQQ[0].c):
    s=fmpq_get_str( NULL, 10, on_row+j )
    S=s
    S=S.replace(' ','')
    if i==0 and j==0:
     r = S
    else:
     r += ','+S
  return r
 
 def raw_str(self):
  '''
  extract all numerators and denominators bypassing fmpq_get_str()
  return result as a space-separated string
  '''
  cdef Py_ssize_t i,j
  cdef Integer num=Integer(0),den=Integer(0)
  cdef mpq_t e
  mpq_init( e ) # no error here, that is the usual way to initialize mpq_t
  cdef fmpq* on_row
  for i in range(self.matQQ[0].r):
   on_row=self.matQQ[0].rows[i]
   for j in range(self.matQQ[0].c):
    fmpq_get_mpq( e, on_row+j )
    mpz_set( num.value, mpq_numref( e ) )
    mpz_set( den.value, mpq_denref( e ) )
    s=raw_entry( num, den )
    if i==0 and j==0:
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
   
   this subroutine translates matrice to fmpz_mat, before doing actual 
    calculation
   So it is a bad idea to call this on big integral matrice
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

def column_to_fmpq_mat( Vector_rational_dense c ):
 ' convert Sage rational vector to fmpq_mat of dimension n*1 '
 cdef fmpq_mat a=fmpq_mat.__new__(fmpq_mat)
 cdef long n=c._degree
 fmpq_mat_init( a.matQQ, n, 1 )
 cdef long i
 for i in range(n):
  fmpq_set_mpq( a.matQQ[0].entries+i, c._entries[i] )
 return a

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

def fmpq_mat_scalar_div(fmpq_mat a, Integer d):
 ' returns a/d as fmpq_mat '
 cdef fmpq_mat b=fmpq_mat.__new__(fmpq_mat)
 fmpq_mat_init( b.matQQ, a.matQQ[0].r, a.matQQ[0].c )
 scalar_div_fmpq_3arg( b, a, d )
 return b

def fmpq_mat_scalar_mul_rational(fmpq_mat a, Rational m):
 ' returns a*m as fmpq_mat '
 cdef fmpq_mat b=fmpq_mat.__new__(fmpq_mat)
 fmpq_mat_init( b.matQQ, a.matQQ[0].r, a.matQQ[0].c )
 cdef fmpq_t q
 fmpq_init( q )
 fmpq_set_mpq( q, m.value )
 cdef long i,size=a.matQQ[0].r * a.matQQ[0].c
 for i in range(size):
  fmpq_mul(  b.matQQ[0].entries+i, a.matQQ[0].entries+i, q )
 fmpq_clear( q )
 return b

def fmpz_mat_mul_by_fmpq_mat(fmpz_mat a, fmpq_mat b):
 '''
 a,b: square matrix
 compute c=a*b
 if c is integer, set a:=c and return 1
 else damage matrice a and return 0
 '''
 cdef fmpq_mat_t C
 cdef long dim=b.matQQ[0].r
 fmpq_mat_init( C, dim, dim )
 fmpq_mat_mul_r_fmpz_mat(C, a.matr, b.matQQ)
 dim=fmpq_mat_get_fmpz_mat(a.matr, C)
 return dim 
