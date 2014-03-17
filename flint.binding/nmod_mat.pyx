# This program is part of RAZIN
# Licence: GNU General Public License (GPL)

'''
nmod_mat_rref and nmod_mat_lu don't do what I want

nmod_mat_lu appears to malfunction for non-prime moduli

Thus this module currently has nothing useful for me
'''

cdef extern from 'flint/nmod_vec.h':
 ctypedef struct nmod_t:
  mp_limb_t n
  mp_limb_t ninv
  mp_bitcnt_t norm

cdef extern from 'flint/nmod_mat.h':
 ctypedef struct nmod_mat_struct:
  mp_limb_t* entries
  long r
  long c
  mp_limb_t** rows
  nmod_t mod
 ctypedef nmod_mat_struct nmod_mat_t[1]
 
 void nmod_mat_init(nmod_mat_t tgt,long rows,long cols,mp_limb_t n)
 void nmod_mat_clear(nmod_mat_t mat)
 long nmod_mat_rref(nmod_mat_t A) 
 # unclear what nmod_mat_rref does. It appears to return spoilt HNF
 long nmod_mat_lu(long *P, nmod_mat_t A, int rank_check)

cdef extern from './nmod_mat_HNF.c':
 void nmod_mat_HNF(nmod_mat_t A)

cdef extern from 'flint/fmpz_mat.h':
 # these two functions undocumented, as of version 2.4.1
 void fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A)
 void fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod)

def AmodB(Integer a, Integer b):
 ' test that mp_limb_t arithmetic works in Cython '
 cdef Integer r=Integer(0)
 cdef mp_limb_t A=mpz_get_ui(a.value)
 cdef mp_limb_t B=mpz_get_ui(b.value)
 cdef mp_limb_t C
 C=A % B
 mpz_set_ui(r.value,C)
 return r

cdef class nmod_mat:
 
 cdef nmod_mat_t matZn
  
 def __init__(self, fmpz_mat z, Integer n):
  '''
   n must be integer between 2 and 2**64-1, else nothing good will happen

   and don't ask me to check n, I won't
  '''
  cdef mp_limb_t nn=mpz_get_ui(n.value)
  nmod_mat_init( self.matZn, z.matr[0].r, z.matr[0].c, nn )
  fmpz_mat_get_nmod_mat(self.matZn, z.matr)

 def __dealloc__(self):
  nmod_mat_clear(self.matZn)
  
 def export_nonnegative_fmpz_mat(self):
  ' lifts oneself to ZZ. returns result as non-negative fmpz_mat '
  cdef Py_ssize_t i,j
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  cdef mp_limb_t* on_row
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  for i in range(self.matZn[0].r):
   on_row=self.matZn[0].rows[i]
   for j in range(self.matZn[0].c):
    fmpz_set_ui( r.matr[0].entries+i*self.matZn[0].c+j, on_row[j] )
  return r

 def export_nonnegative_fmpz_mat_upper(self):
  ' like export_nonnegative_fmpz_mat, but only outputs diagonal and upper part '
  cdef Py_ssize_t i,j
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  cdef mp_limb_t* on_row
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  for i in range(self.matZn[0].r):
   on_row=self.matZn[0].rows[i]
   for j in range(i,self.matZn[0].c):
    fmpz_set_ui( r.matr[0].entries+i*self.matZn[0].c+j, on_row[j] )
  return r

 def export_balanced_fmpz_mat(self):
  '''
  lifts oneself to ZZ. returns result as fmpz_mat with possibly negative entries
  
  this function currently untested
  '''
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  fmpz_mat_set_nmod_mat( r.matr, self.matZn )
  return r

 def rref(self):
  '''
  Bring oneself into row-echelon form 
  
  return rank
  
  I don't know what this procedure does, but it certainly does not count HNF
   because it does not preserve deteminant
  '''
  cdef Integer r=Integer(0)
  cdef long rank
  rank=nmod_mat_rref(self.matZn)
  mpz_set_si(r.value,rank)
  return r

 def __repr__(self):
   return 'nmod_mat(%i, %i, [%s])' % (self.matZn[0].r, self.matZn[0].c,
            (', '.join(map(str, self.entries()))))

 def entries(self):
  cdef Py_ssize_t i,j
  cdef Integer t
  cdef mp_limb_t* on_row
  r,t=[],Integer(0)
  for i in range(self.matZn[0].r):
   on_row=self.matZn[0].rows[i]
   for j in range(self.matZn[0].c):
    mpz_set_ui( t.value, on_row[j] )
    r.append( int(t) )
  return r

def fmpz_mat_hermite_form(fmpz_mat A,Integer M):
 '''
 A: square non-singular over Z matrice
 
 M: in range 2..2**64-1, det(A) divides M
 
 return Pernet/Stein/Sage HNF of A over ring of integers as fmpz_mat

 Before asking me questions, count HNF of 
  [  5  5 ]
  [ -7 -8 ]
 and 
  [   66551115 -1211111295]
  [  -92996400  1692368196] 
 using modular technique
 '''
 a=nmod_mat(A,M)
 #sig_on() # added temporarily to shorten assert failed dump
 nmod_mat_HNF(a.matZn)
 #sig_off()
 return a.export_nonnegative_fmpz_mat_upper()

def fmpz_mat_hermite_form____not_working(fmpz_mat A,Integer M):
 '''
 A: n*m matrice with m>=n, whose 1st n columns form a non-singular 
  matrice B
 M: integer in range 2..2**64-1, det(B) divides M
 return Pernet/Stein/Sage HNF of A over ring of integers as fmpz_mat
 
 ups, for non-prime M this subroutine does not work
 '''
 a=nmod_mat(A,M)
 print 'taken modulo',M,':',a
 cdef long* p=<long*>malloc(a.matZn[0].r * sizeof(long))
 nmod_mat_lu(p,a.matZn,<int>0)
 '''
 for matrice 5,5,3,2 modulo 10 nmod_mat_lu returns 5,5,3,7 
 
 but equality 
 1 0   5 5   5 5
 3,1 * 0 7 = 3 2
 does not hold
 
 However, for same matrice 5,5,-7,-8 modulo 13 result is correct:
 
 5 5   1 0   5  5
 6 5 = 9 1 * 0 12
 '''
 print 'hermite form:',a
 print 'p',p[0],' ',p[1]
 cdef fmpz_mat b=a.export_nonnegative_fmpz_mat_upper()
 free(p)
 return b
