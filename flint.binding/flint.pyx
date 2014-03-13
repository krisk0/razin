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
from sage.rings.rational cimport Rational

from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense

include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"

# numbers imported
cdef extern from 'flint/fmpz.h':
 ctypedef long fmpz_t[1]
 void fmpz_set_mpz(fmpz_t tgt, mpz_t sou)
 void fmpz_init(fmpz_t x)
 void fmpz_get_mpz(mpz_t tgt, fmpz_t sou)
 void fmpz_clear(fmpz_t f)

cdef extern from 'flint/fmpq.h':
 ctypedef struct fmpq:
  long num
  long den
 ctypedef fmpq fmpq_t[1]
 char* fmpq_get_str(char* rez, int b, const fmpq_t x)
 
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

# Python types and methods
include "fmpz_mat.pxi"
include "fmpq_mat.pxi"
