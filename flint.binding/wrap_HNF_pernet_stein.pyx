# This program is part of RAZIN

# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# cython definitions required to wrap fmpz_mat_hnf_pernet_stein() 

cdef extern from 'C/fmpz_mat/hnf_pernet_stein.c':
 void fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A)

def hnf_pernet_stein(fmpz_mat a):
 cdef fmpz_mat_t h
 fmpz_mat_init(h, a.matr.r, a.matr.c)
 fmpz_mat_hnf_pernet_stein(h, a.matr)
 return wrap_fmpz_mat(h)

cdef extern from 'C/fmpz_mat/hnf_modular.c':
 void fmpz_mat_hnf_modular(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D)

cdef extern from 'C/fmpz_mat/hnf_xgcd.c':
 void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A)

cdef extern from 'C/fmpz_mat/hnf_classical.c':
 void fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A)
 
cdef extern from 'C/fmpz_mat/hnf_minors.c':
 void fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A)
