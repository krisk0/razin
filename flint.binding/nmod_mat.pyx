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
 
 long nmod_mat_rref(nmod_mat_t A)
 void nmod_mat_clear(nmod_mat_t mat)
 void nmod_mat_init(nmod_mat_t tgt,long rows,long cols,mp_limb_t n)

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
  cdef Py_ssize_t size,i
  cdef fmpz_t t
  fmpz_init( t )
  nmod_mat_init( self.matZn, z.matr[0].r, z.matr[0].c, nn )
  size = z.matr[0].r * z.matr[0].c
  for i in range(size):
   self.matZn[0].entries[i]=fmpz_mod_ui( t, z.matr[0].entries+i, nn )
  fmpz_clear( t )

 def __dealloc__(self):
  nmod_mat_clear(self.matZn)
  
 def export_fmpz_mat(self):
  ' lifts oneself to ZZ. returns result as fmpz_mat '
  cdef Py_ssize_t size,i
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  size = self.matZn[0].r * self.matZn[0].c
  for i in range(size):
   fmpz_set_ui( r.matr[0].entries+i, self.matZn[0].entries[i] )
  return r
