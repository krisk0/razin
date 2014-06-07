# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

# AKS code put into public domain by 袁轶君 (Yijun Yuan), then rewritten completely
#  by me
cdef extern from 'C/fmpz/AKS_trunc.c':
 int AKS_ui(mp_limb_t n)
 int AKS(fmpz_t n)

def n_is_prime_AKS(n):
 'n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64'
 return AKS_ui(<mp_limb_t>n)

cdef extern from 'C/ulong_extras/big_primes_iterator.c':
 int n_is_big_sureprime(mp_limb_t n)

cdef extern from 'C/ulong_extras/gcd_odd_.c':     # this only compiles on 
 mp_limb_t n_gcd_odd_even(mp_limb_t x,mp_limb_t y)#  AMD/Intel
 mp_limb_t n_gcd_odd_odd(mp_limb_t x,mp_limb_t y)

cdef extern from 'C/ulong_extras/inv_mod_pk.c':
 # according to FLINT documentation mp_limb_t is same size and unsigned like 
 #  ulong. Can't write ulong on Cython, replacing with mp_limb_t
 mp_limb_t inv_mod_pk(mp_limb_t a,mp_limb_t p,mp_limb_t k,mp_limb_t p_deg_k,
  mp_limb_t p_deg_k_norm,mp_limb_t p_deg_k_inv)

cdef extern from 'flint/ulong_extras.h':
 int n_is_probabprime_BPSW(mp_limb_t n)
 mp_limb_t n_nextprime(mp_limb_t n, int proved)
 mp_limb_t n_preinvert_limb(mp_limb_t n)
 
def is_big_sureprime(n):
 return n_is_big_sureprime(<mp_limb_t>n)

def is_probabprime_BPSW(n):
 return n_is_probabprime_BPSW(<mp_limb_t>n)

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

def prev_prime(a):
 '''
  return previous prime for a >= MIN_n_primes_rev
  else return next prime for a-1
 '''
 cdef mp_limb_t n=a, p
 if n&1:
  n += 1
 cdef n_primes_rev_t i
 p=n_primes_rev_init(i,n)
 n_primes_rev_clear(i)
 if p>1:
  return p
 return n_nextprime( (<mp_limb_t>a)-1, 0)

def count_primes_in_range(a,b):
 '''
  MIN_n_primes_rev <= a 
  (a < b < 2**64, b prime or even) or b=0
  
  return count of prime numbers in range a..b' where b'=b or 0xFFFFFFFFFFFFFFC5
 '''
 cdef mp_limb_t aa=a, bb=b, p, r=0
 cdef n_primes_rev_t i
 p=n_primes_rev_init(i,b)
 while 1:
  r += 1
  p=n_primes_rev_next(i)
  if p<aa:
   break
 n_primes_rev_clear(i)
 return int(r)

def primes_in_range_2010(a,b):
 '''
 find all primes in range a..b using FLINT n_nextprime()
 
 see notes in count_primes_in_range_2010() below on range of a and b
 '''
 cdef mp_limb_t aa=a, bb=b, p
 if aa & 1:
  aa -= 1
 p=n_nextprime(aa,0)  # even though 2nd parameter proved is 0, error is unlikely
 r=dict()
 while 1:
  if p>b:
   break
  r[ int(p) ] = 2
  p=n_nextprime(p,0)
 return r

def count_primes_in_range_2010(a,b):
 '''
 return count of prime numbers in range a..b, valid range for a and b is not
  known
 
 according to FLINT documentation, result is guaranteed for 2<a<b<10**16
 
 b >= 0xFFFFFFFFFFFFFFC5 results in exception
 '''
 cdef mp_limb_t aa=a, bb=b, p, r=0
 if aa & 1:
  aa -= 1
 p=n_nextprime(aa,0)
 while 1:
  if p>b:
   break
  r += 1
  p=n_nextprime(p,0)
 return r

def inverse_mod_pk(a,p,k,p_deg_k,p_deg_k_nrm,p_deg_k_inv):
 return inv_mod_pk(<mp_limb_t>a,<mp_limb_t>p,<mp_limb_t>k,
  <mp_limb_t>p_deg_k,<mp_limb_t>p_deg_k_nrm,<mp_limb_t>p_deg_k_inv)

def inverse_mod_pk_3arg(a,p,k):
 cdef mp_limb_t p_deg_k = <mp_limb_t>( p**k )
 cdef p_deg_k_nrm=p_deg_k
 while p_deg_k_nrm<0x8000000000000000:
  p_deg_k_nrm <<= 1
 cdef p_deg_k_inv=n_preinvert_limb(p_deg_k_nrm)
 return inv_mod_pk(<mp_limb_t>a,<mp_limb_t>p,<mp_limb_t>k,
  p_deg_k,p_deg_k_nrm,p_deg_k_inv)

def n_preinvert_limb_wr( n ):
 return n_preinvert_limb(<mp_limb_t>n)
