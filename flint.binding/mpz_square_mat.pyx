# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

cdef extern from 'C/mpz_square_mat/mul_vec_mat_modulo.c':
 void mpz_square_mat_mul_vec_mat_modulo(mpz_ptr t,mpz_ptr v,
  const mpz_square_mat_t A,mpz_t m)
