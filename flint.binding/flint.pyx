# -*- coding: utf-8

# This program is part of RAZIN

###########################################################################

 #original code in .../sage/libs/flint
#       Copyright (C) 2013 Fredrik Johansson <fredrik.johansson@gmail.com>
 #modified by Денис Крыськов to access some fmpX_mat_XXX() subroutines 

#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

# I am adding new functions. Only functions I personally need to compute HNF
# Maybe To.do file will give a clue, on what and why I am doing

from sage.rings.integer cimport Integer
#from sage.rings.integer_ring import ZZ
from sage.rings.rational cimport Rational
from sage.rings.rational_field import QQ
from sage.matrix.constructor import Matrix

from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.modules.vector_rational_dense cimport Vector_rational_dense
from sage.modules.free_module_element import vector

from libc.stdlib cimport malloc, free
include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"

# numbers imported
cdef extern from 'gmp.h':
 ctypedef unsigned long mp_limb_t
 ctypedef unsigned long mp_bitcnt_t

ctypedef __mpz_struct* mpz_ptr

cdef extern from 'flint/fmpz.h':
 ctypedef long fmpz_t[1]
 void fmpz_set_mpz(fmpz_t tgt, mpz_t sou)
 void fmpz_init(fmpz_t x)
 void fmpz_get_mpz(mpz_t tgt, fmpz_t sou)
 void fmpz_clear(fmpz_t f)
 mp_limb_t fmpz_mod_ui(fmpz_t f, const fmpz_t g, mp_limb_t x)
 void fmpz_set_ui(fmpz_t tgt, mp_limb_t val)
 void fmpz_lcm(fmpz_t f , const fmpz_t g , const fmpz_t h)
 void fmpz_one(fmpz_t tgt)
 
cdef extern from 'flint/fmpq.h':
 ctypedef struct fmpq:
  long num
  long den
 ctypedef fmpq fmpq_t[1]
 void fmpq_init( fmpq_t x )
 void fmpq_clear( fmpq_t x )
 char* fmpq_get_str(char* rez, int b, const fmpq_t x)
 void fmpq_set(fmpq_t tgt, const fmpq_t src)
 void fmpq_set_fmpz_frac(fmpq_t tgt, const fmpz_t p, const fmpz_t q)
 void fmpq_set_mpq(fmpq_t tgt, const mpq_t src)
 void fmpq_get_mpq(mpq_t tgt, const fmpq_t src)
 void fmpq_div_fmpz(fmpq_t tgt, const fmpq_t sou, const fmpz_t x)
 void fmpq_mul(fmpq_t tgt, const fmpq_t a, const fmpq_t b)

# matrix imported
cdef extern from 'flint/fmpz_mat.h':
 ctypedef struct fmpz_mat_struct:
  long* entries
  long r
  long c
  long** rows
 ctypedef fmpz_mat_struct fmpz_mat_t[1]

cdef extern from 'flint/fmpq_mat.h':
 ctypedef struct fmpq_mat_struct:
  fmpq* entries
  long r
  long c
  fmpq** rows
 ctypedef fmpq_mat_struct fmpq_mat_t[1]

# Python-visible types and methods
include "fmpz_mat.pyx"
include "fmpq_mat.pyx"
include "nmod_mat.pyx"
include "tmod_mat.pyx"
