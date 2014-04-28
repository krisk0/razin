# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

# AKS code put into public domain by 袁轶君 (Yijun Yuan)

cdef extern from 'C/fmpz/AKS_trunc.c':
 int AKS_ui(mp_limb_t n)
 int AKS(fmpz_t n)

def n_is_prime_AKS(Integer n):
 'n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64'
 cdef mp_limb_t nn=mpz_get_ui(n.value)
 return AKS_ui(nn)
