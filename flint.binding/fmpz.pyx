# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

# AKS code put into public domain by 袁轶君 (Yijun Yuan), then rewritten completely
#  by me
cdef extern from 'C/fmpz/AKS_trunc.c':
 int AKS_ui(mp_limb_t n)
 int AKS(fmpz_t n)

def n_is_prime_AKS(Integer n):
 'n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64'
 cdef mp_limb_t nn=mpz_get_ui(n.value)
 return AKS_ui(nn)

cdef extern from 'C/ulong_extras/big_primes_iterator.c':
 int n_is_big_sureprime(mp_limb_t n)

cdef extern from 'flint/ulong_extras.h':
 int n_is_probabprime_BPSW(mp_limb_t n)

def is_big_sureprime(n):
 cdef mp_limb_t nn=n
 return n_is_big_sureprime(nn)

def is_probabprime_BPSW(n):
 cdef mp_limb_t nn=n
 return n_is_probabprime_BPSW(nn)

cdef extern from 'C/ulong_extras/big_primes_iterator.c':
 ctypedef struct n_primes_rev_struct:
  long index
  long allocated_size
  long last_output_mod_30
  mp_limb_t* numbers
 ctypedef n_primes_rev_struct n_primes_rev_t[1]
 mp_limb_t n_primes_rev_init(n_primes_rev_t i, mp_limb_t start)
 void n_primes_rev_clear(n_primes_rev_t i)
 mp_limb_t n_primes_rev_show_again(n_primes_rev_t i)
 mp_limb_t n_primes_rev_next(n_primes_rev_t i)
 mp_limb_t n_primes_rev_reset(n_primes_rev_t i)

def primes_in_range(a,b):
 '''
  MIN_n_primes_rev <= a 
  (a < b < 2**64, b prime or even) or b=0
  
  return prime numbers in range a..b' where b'=b or 0xFFFFFFFFFFFFFFC5 
  
  result in form of Python dict, prime->1
 '''
 cdef mp_limb_t aa=a, bb=b, p
 cdef n_primes_rev_t i
 p=n_primes_rev_init(i,b)
 r=dict()
 while 1:
  r[ int(p) ]=1
  p=n_primes_rev_next(i)
  if p<aa:
   break
 n_primes_rev_clear(i)
 return r
