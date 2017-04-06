# This program is part of RAZIN

# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# cython definitions required to wrap HNF subroutines 

cdef extern from 'flint/flint.h':
 void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)

def AlexBest_hnf(fmpz_mat a):
 cdef fmpz_mat_t h
 fmpz_mat_init(h, a.matr.r, a.matr.c)
 fmpz_mat_hnf(h, a.matr)
 return wrap_fmpz_mat(h)
