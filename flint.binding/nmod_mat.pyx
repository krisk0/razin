# This program is part of RAZIN
# Licence: GNU General Public License (GPL)

cdef extern from "gmp.h":
 ctypedef unsigned long mp_limb_t
 ctypedef unsigned long mp_bitcnt_t

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
