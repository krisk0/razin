# This program is part of RAZIN
# Licence: GNU General Public License (GPL)

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
  cdef Py_ssize_t size,i
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  size = self.matZn[0].r * self.matZn[0].c
  for i in range(size):
   fmpz_set_ui( r.matr[0].entries+i, self.matZn[0].entries[i] )
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

 def hermite_form(self):
  '''
  Bring oneself into row-echelon form aka Pernet/Stein/Sage Hermite Normal Form
  
  return rank
  '''
  cdef Integer r=Integer(0)
  cdef long rank
  rank=nmod_mat_rref(self.matZn)
  mpz_set_si(r.value,rank)
  return r

def fmpz_mat_hermite_form(fmpz_mat A,Integer M):
 '''
 A: n*m matrice with m>=n, whose 1st n columns form a non-singular 
  matrice B
 M: integer in range 2..2**64-1, such that det(B) divides M
 
 return Pernet/Stein/Sage HNF of A over ring of integers as fmpz_mat
 '''
 a=nmod_mat(A,M)
 rank=a.hermite_form()
 cdef fmpz_mat b=a.export_nonnegative_fmpz_mat()
 if rank<a.matZn[0].r:
  ' if rank over residue ring is smaller, then put M into lower-right corner '
  fmpz_set_mpz( b.matr[0].entries + (a.matZn[0].r * a.matZn[0].c - 1), M.value )
 return b
