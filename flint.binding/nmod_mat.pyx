# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

cdef extern from 'flint/nmod_mat.h':
 void nmod_mat_init(nmod_mat_t tgt,long rows,long cols,mp_limb_t n)
 void nmod_mat_clear(nmod_mat_t mat)
 long nmod_mat_rref(nmod_mat_t A) 
 # unclear what nmod_mat_rref does. It appears to return spoilt HNF
 long nmod_mat_lu(long *P, nmod_mat_t A, int rank_check)
 mp_limb_t nmod_mat_det(nmod_mat_t A)

cdef extern from 'C/tmod_mat/tmod_mat.c':
 mp_limb_t fmpz_to_t(const fmpz_t f)

cdef extern from 'C/nmod_mat/det_mod_pk.c':
 mp_limb_t nmod_mat_det_mod_pk(nmod_mat_t M,mp_limb_t p,mp_limb_t k,
  mp_limb_t p_deg_k,mp_limb_t* scratch)
 mp_limb_t det_mod_pk_fix_SE_corner(mp_limb_t* I,nmod_mat_t M,long* negate_det,
  mp_limb_t p,mp_limb_t k,mp_limb_t p_deg_k)

cdef extern from 'C/nmod_mat/init_square_2arg.c':
 void nmod_mat_init_square_2arg(nmod_mat_t mat, long dim)

cdef extern from 'flint/fmpz_mat.h':
 # these two functions undocumented, as of version 2.4.1
 void fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A)
 void fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod)
 void fmpz_mat_zero(fmpz_mat_t A)

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
  nmod_mat_init( self.matZn, z.matr[0].r, z.matr[0].c, nn )
  fmpz_mat_get_nmod_mat(self.matZn, z.matr)

 def __dealloc__(self):
  nmod_mat_clear(self.matZn)
  
 def export_nonnegative_fmpz_mat(self):
  ' lifts oneself to ZZ. returns result as non-negative fmpz_mat '
  cdef Py_ssize_t i,j
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  cdef mp_limb_t* on_row
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  for i in range(self.matZn[0].r):
   on_row=self.matZn[0].rows[i]
   for j in range(self.matZn[0].c):
    fmpz_set_ui( r.matr[0].entries+i*self.matZn[0].c+j, on_row[j] )
  return r

 def export_nonnegative_fmpz_mat_upper(self):
  ' like export_nonnegative_fmpz_mat, but only outputs diagonal and upper part '
  return nmod_mat_export_nonnegative_fmpz_mat_upper( self.matZn )
  
 def export_balanced_fmpz_mat(self):
  '''
  lifts oneself to ZZ. returns result as fmpz_mat with possibly negative entries
  
  this function currently untested
  '''
  cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
  fmpz_mat_init( r.matr, self.matZn[0].r, self.matZn[0].c )
  fmpz_mat_set_nmod_mat( r.matr, self.matZn )
  return r

 def rref(self):
  '''
  Bring oneself into row-echelon form 
  
  return rank
  
  I don't know what this procedure does, but it certainly does not count HNF
   because it does not preserve deteminant
  '''
  cdef Integer r=Integer(0)
  cdef long rank
  rank=nmod_mat_rref(self.matZn)
  mpz_set_si(r.value,rank)
  return r

 def __repr__(self):
   return 'nmod_mat(%i, %i, [%s])' % (self.matZn[0].r, self.matZn[0].c,
            (', '.join(map(str, self.entries()))))

 def entries(self):
  cdef Py_ssize_t i,j
  cdef Integer t
  cdef mp_limb_t* on_row
  r,t=[],Integer(0)
  for i in range(self.matZn[0].r):
   on_row=self.matZn[0].rows[i]
   for j in range(self.matZn[0].c):
    mpz_set_ui( t.value, on_row[j] )
    r.append( int(t) )
  return r
  
 def determinant(self):
  cdef Integer r=Integer(0)
  cdef mp_limb_t d=nmod_mat_det(self.matZn)
  mpz_set_ui( r.value, d )
  return r

def det_modulo_prime(Matrix_integer_dense sage, Integer p):
 # TODO: to speed-up, directly convert Sage matrice to nmod_mat
 cdef fmpz_mat a=fmpz_mat( sage )
 return nmod_mat( a, p ).determinant()

cdef nmod_mat_export_nonnegative_fmpz_mat_upper( nmod_mat_t s ):
 ' export upper part of s as fmpz_mat '
 cdef Py_ssize_t i,j
 cdef Py_ssize_t cc=s.c, rc=s.r
 cdef fmpz_mat r=fmpz_mat.__new__(fmpz_mat)
 fmpz_mat_init( r.matr, rc, cc )
 cdef mp_limb_t* on_s_row
 cdef long* on_r_row
 for i in range(rc):
  on_s_row=s.rows[i]
  on_r_row=r.matr[0].rows[i]
  for j in range(i,cc):
   fmpz_set_ui( on_r_row+j, on_s_row[j] )
 return r

cdef fmpz_mat_set_nmod_mat_upper(fmpz_mat_t tgt, nmod_mat_t sou):
 '''
 cc <= rc, tgt and sou has same dimenisions. 
 Zero tgt, make it virgin, put upper part of sou into tgt 
 '''
 cdef Py_ssize_t i,j
 cdef Py_ssize_t cc=sou.c, rc=sou.r
 fmpz_mat_zero(tgt)
 cdef mp_limb_t* on_s_row
 cdef long* on_r_row=tgt.entries
 for i in range(cc):
  on_s_row=sou.rows[i]
  tgt.rows[i]=on_r_row
  for j in range(i,cc):
   fmpz_set_ui( on_r_row+j, on_s_row[j] )
  on_r_row += cc

def fmpz_mat_hermite_form(fmpz_mat A, Integer M):
 '''
 A: square non-singular over Z matrice
 
 M: in range 2..2**64-1, abs(det(A)) divides M
 
 return Pernet/Stein/Sage HNF of A over ring of integers as fmpz_mat

 Before asking me questions, count HNF of 
  [  5  5 ]
  [ -7 -8 ]
 modulo 5 and 
  [   66551115 -1211111295]
  [  -92996400  1692368196] 
 modulo 540
 '''
 a=nmod_mat(A,M)
 #sig_on() # added temporarily to shorten assert failed dump
 nmod_mat_HNF(a.matZn)
 #sig_off()
 return a.export_nonnegative_fmpz_mat_upper()

cdef export_triU(fmpz_mat_t R, nmod_mat_t S):
 '''
 on entry R un-initilized, S.r >= S.c
 on exit R: square upper-trianglular, formed from upper part of S (lower part
  ignored)
 '''
 cdef long dim=S.c,i,j
 fmpz_mat_init( R, dim, dim )
 cdef long* on_tgt
 cdef mp_limb_t* on_src
 for i in range(dim):
  on_src = S.rows[i]+i
  on_tgt = R.rows[i]+i
  for j in range(dim-i):
   fmpz_set_ui( on_tgt, on_src[0] )
   on_tgt += 1
   on_src += 1

def fmpz_mat_hermite_form_nonsquare(fmpz_mat A, Integer M):
 '''
 A rows count >= A columns count=c,    1 < M < 2**64
 A non-singular over Z, M divides det( H' ) where H' = upper c rows of HNF(A)
 call nmod_mat_HNF_nonsquare()
 return H' as fmpz_mat
 '''
 a=nmod_mat(A, M)
 nmod_mat_HNF_nonsquare( a.matZn )
 cdef fmpz_mat_t r
 export_triU(r, a.matZn)
 return wrap_fmpz_mat( r )

def plu_modulo_prime(fmpz_mat A,Integer M):
 '''
 call nmod_mat_lu on A taken modulo M, get P and jammed LU 
 
 set B=P*A
 
 return B,LU as fmpz_mat
 '''
 a=nmod_mat(A,M)
 cdef long* p=<long*>malloc(a.matZn[0].r * sizeof(long))
 nmod_mat_lu(p,a.matZn,<int>0)
 '''
 for matrice 5,5,3,2 modulo 10 nmod_mat_lu returns 5,5,3,7 
 
 but equality 
 1 0   5 5   5 5
 3,1 * 0 7 = 3 2
 does not hold
 
 However, for same matrice 5,5,-7,-8 modulo 13 result is correct:
 
 5 5   1 0   5  5
 6 5 = 9 1 * 0 12
 '''
 cdef fmpz_mat LU=a.export_nonnegative_fmpz_mat()
 cdef fmpz_mat B=fmpz_mat_permute( p, A )
 free(p)
 return B,LU

def test_invert_4x4_corner(fmpz_mat A,p,k):
 '''
 call det_mod_pk_fix_SE_corner
 if result less than 3, return None
 else return tuple (negate_det,determinant,Mi,M)
 '''
 cdef mp_limb_t pp=p
 cdef mp_limb_t kk=k
 cdef mp_limb_t p_deg_k = <mp_limb_t>( p**k )
 cdef Integer p_deg_k_nrm=Integer(0)
 mpz_set_ui( p_deg_k_nrm.value, p_deg_k )
 while p_deg_k_nrm<0x8000000000000000:
  p_deg_k_nrm <<= 1
 M=nmod_mat(A,p_deg_k_nrm)
 cdef nmod_mat Mi=nmod_mat.__new__( nmod_mat )
 nmod_mat_init( Mi.matZn, 4, 4, p_deg_k_nrm )
 cdef long negate_det[1]
 negate_det[0]=0
 cdef mp_limb_t r=det_mod_pk_fix_SE_corner(Mi.matZn.entries,M.matZn,negate_det,
  pp,kk,p_deg_k)
 if r<3:
  return
 r = (r-2) % p_deg_k
 return negate_det[0],r,Mi.export_nonnegative_fmpz_mat().export_sage(),\
  M.export_nonnegative_fmpz_mat().export_sage()

cdef void nmod_mat_mod_t_half(nmod_mat_t tgt, fmpz_mat_t sou):
 '''
 compute tgt := sou modulo 2**63
 
 This subroutine is for amd64
 '''
 cdef long rc=sou.r, cc=sou.c, i, j
 nmod_mat_init( tgt, rc, 0, 0x8000000000000000 )
 cdef mp_limb_t* e = <mp_limb_t*>flint_malloc( rc * cc * sizeof(mp_limb_t) )
 tgt.rows = <mp_limb_t**>flint_malloc( rc * sizeof(mp_limb_t*) )
 cdef long* s
 tgt.entries = e
 tgt.c=cc
 for i in range(rc):
  tgt.rows[i] = e
  s=sou.rows[i]
  for j in range(cc):
   e[j] = fmpz_to_t( s+j ) & 0x7FFFFFFFFFFFFFFF
  e += cc

def nmod_mat_set_fmpz_mat_mod_thalf(fmpz_mat A):
 cdef nmod_mat r=nmod_mat.__new__( nmod_mat )
 nmod_mat_mod_t_half( r.matZn, A.matr )
 return r
