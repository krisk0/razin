# -*- coding: utf-8
# This program is part of RAZIN
# Licence: GNU General Public License (GPL)
# Copyright Денис Крыськов 2014

cdef extern from 'flint/flint.h':
 ctypedef struct flint_rand_s:
  gmp_randstate_t gmp_state
  int gmp_init
  mp_limb_t __randval
  mp_limb_t __randval2
 ctypedef flint_rand_s flint_rand_t[1]
 void flint_randinit(flint_rand_t state)
 void flint_randclear(flint_rand_t state)
 ctypedef unsigned long mp_bitcnt_t
 void fmpz_randbits(fmpz_t tgt, flint_rand_t state, mp_bitcnt_t bits)
 void fmpz_randtest(fmpz_t tgt, flint_rand_t state, mp_bitcnt_t bits)

cdef extern from './C/tmod_mat/tmod_mat.c':
 ctypedef struct tmod_mat_struct:
  long r
  long c
  mp_limb_t** rows
  mp_limb_t* entries
 ctypedef tmod_mat_struct tmod_mat_t[1]
 mp_limb_t fmpz_to_t(const fmpz_t f)
 void tmod_mat_init(tmod_mat_t m, long rows, long cols)
 void tmod_mat_init_fast(tmod_mat_t mat, long rows, long cols)
 void tmod_mat_clear(tmod_mat_t m)
 long tmod_mat_PLU_mod_machine_word(long* PR, tmod_mat_t S)
 void tmod_mat_solver_3arg( tmod_mat_t R, long* PD, const tmod_mat_t LU )

cdef mp_limb_t test_fmpz_to_t(fmpz_t n,flint_rand_t S,mp_bitcnt_t bits,
  long upto):
 ' return 0 iff test passes '
 if upto:
  fmpz_randtest(n,S,bits)
 else:
  fmpz_randbits(n,S,bits)
 cdef Integer m=Integer(0)
 fmpz_get_mpz( m.value, n )
 cdef mp_limb_t k0,k1
 k0=<mp_limb_t>( m % 2**64 )
 k1=fmpz_to_t(n)
 if k0==k1:
  return 0
 print 'test failed for n=%s, k0=%X k1=%X' % (m,k0,k1)
 return 1

def check_fmpz_to_t():
 '''
 test fmpz_get_junior_limb:
  generate a random number as described below
  extract junior limb with fmpz_get_junior_limb()
  convert to Integer, extract junior limb by taking it 2**64
  results should match
  
 return 1 iff test passed
  
 Random numbers choice: 
  1000 up to 32 bits
  1000 up to 62 bits
  1000 of bit-length 61
  1000 of bit-length 62
  1000 of bit-length 63
  1000 of bit-length 64
  100 of bit-length 65
  100 of bit-length 66
  100 of bit-length 127
  100 of bit-length 128
  100 of bit-length 129
 '''
 cdef fmpz_t t
 fmpz_init( t )
 cdef flint_rand_t S
 flint_randinit(S)
 cdef long i
 cdef mp_bitcnt_t j
 for i in range(1000):
  for j in 32,62:
   if test_fmpz_to_t(t,S,j,1):
    return 0
 for i in range(1000):
  for j in range(61,65):
   if test_fmpz_to_t(t,S,j,0):
    return 0
 for i in range(100):
  for j in 65,66,127,128,129:
   if test_fmpz_to_t(t,S,j,0):
    return 0
 flint_randclear(S)
 fmpz_clear( t )
 return 1

def take_mod_2_64(Integer x):
 '''
 assist in testing fmpz_to_t() subroutine from Python
 
 convert argument to fmpz, then call fmpz_to_t, return its result
 '''
 cdef fmpz_t d
 fmpz_init( d )
 fmpz_set_mpz(d, x.value)
 cdef mp_limb_t r=fmpz_to_t(d)
 fmpz_clear( d )
 return r

cdef tmod_mat_set_fmpz_mat( tmod_mat_t tgt, fmpz_mat_t matr ):
 tmod_mat_init_fast( tgt, matr[0].r, matr[0].c )
 cdef long i,j,c=matr[0].c,k=0
 cdef long* on_row
 for i in range(matr[0].r):
  on_row=matr[0].rows[i]
  for j in range(c):
   tgt.entries[k]=fmpz_to_t( on_row+j )
   k += 1

def tmod_mat_permute_fmpz_mat( agnostic_array p, fmpz_mat src ):
 '''
 r=src converted t_mod_mat
 apply p to rows of r
 return r
 '''
 cdef tmod_mat_t tgt
 tmod_mat_init( tgt, src.matr[0].r, src.matr[0].c )
 cdef long i,j,c=tgt.c
 cdef long*      on_src_row
 cdef mp_limb_t* on_tgt_row
 cdef long* P=<long*>p.array
 for i in range( tgt.r ):
  on_src_row = src.matr[0].rows[ P[i] ]
  on_tgt_row =         tgt.rows[   i  ]
  for j in range(c):
   on_tgt_row[j]=fmpz_to_t( on_src_row+j )
 return wrap_tmod_mat( tgt )

cdef class tmod_mat_single:
 
 cdef tmod_mat_t matT
  
 def __init__(self, fmpz_mat z):
  tmod_mat_set_fmpz_mat( self.matT, z.matr )

 def __dealloc__(self):
  tmod_mat_clear(self.matT)

 def export_sage(self):
  ' export self as sage Matrix_integer_dense '
  cdef Matrix_integer_dense r=Matrix( self.matT[0].r, self.matT[0].c )
  cdef Py_ssize_t i,j,k=0
  cdef mp_limb_t* on_row
  for i in range(self.matT[0].r):
   on_row=self.matT[0].rows[i]
   for j in range(self.matT[0].c):
    mpz_set_ui( <mpz_ptr>(r._entries+k), on_row[j] )
    k += 1
  return r
 
def export_L_sage( tmod_mat_single me ):
 '''
 un-compress L part of matrice, return result as Sage Matrix_integer_dense,
  square
 
 count of rows of me should be equal or 1 more than count of columns
 '''
 cdef Py_ssize_t i,j, r=me.matT[0].r
 cdef Matrix_integer_dense L=Matrix( r, r )
 cdef mp_limb_t* on_s_row
 cdef mpz_t*     on_t_row
 mpz_set_ui( <mpz_ptr>(L._entries+0), 1 )
 for i in range( 1, r ):
  on_s_row = me.matT[0].rows[i]
  on_t_row = L._entries + i * r
  mpz_set_ui( <mpz_ptr>on_t_row[i], 1 )
  for j in range(i):
   mpz_set_ui( <mpz_ptr>on_t_row[j], on_s_row[j] )
 return L
  
def export_U_sage( tmod_mat_single me ):
 '''
 un-compress U part of matrice, return result as Sage Matrix_integer_dense

 count of rows or me should be not less than count of columns
 '''
 cdef Py_ssize_t i,j,c=me.matT[0].c
 cdef Matrix_integer_dense r=Matrix( me.matT[0].r, c )
 cdef mp_limb_t* on_s_row
 cdef mpz_t*     on_t_row
 for i in range( c ):
  on_s_row=me.matT[0].rows[i]
  on_t_row=r._entries+i*c
  for j in range( i, c):
   mpz_set_ui( <mpz_ptr>(on_t_row+j), on_s_row[j] )
 return r

cdef wrap_tmod_mat(tmod_mat_t a):
 '''
 shallow constructor: wrap existing tmod_mat_t struct as tmod_mat_single Python
  object
 '''
 cdef tmod_mat_single A=tmod_mat_single.__new__( tmod_mat_single )
 A.matT.entries=a.entries
 A.matT.r=a.r
 A.matT.c=a.c
 A.matT.rows=a.rows
 return A

def tmod_mat_PLU(fmpz_mat s):
 '''
 s:matrice of dim n*(n-1), n-1 >= 1, whose HNF modulo 2 is diagonal of ones
  and one zero line
 
 on error return None,None
 
 on success return pair PR,LU where 
  PR is agnostic_array of long representing P and R:
   P: lower line of permutation that need to be applied to rows of s
   R: diagonal of U inverted
  LU: single matrice representing two matrices L and U, as 
   FLINT nmod_mat_lu_classical() does
  
 P * old s = L * U 
 
  1  -6         1  -6               -4 -15
  4 -21   ->    4   3     (inv->)    1   2
  7 -48         7  -2                2 1/3
 
 '''
 cdef tmod_mat_t m
 # delay Python object creation for m
 tmod_mat_set_fmpz_mat( m, s.matr )
 cdef void* PR=malloc( (m.r+m.c) * sizeof(long) )
 if( tmod_mat_PLU_mod_machine_word( <long*>PR, m ) ):
  return wrap_agnostic_array( PR ),wrap_tmod_mat( m )
 # not unimodular matrice. Cleanup and return None
 tmod_mat_clear(m)
 return None,None

def tmod_mat_solver(agnostic_array PD, tmod_mat_single LU):
 '''
 PD, LU: jammed PR, jammed LU as returned by tmod_mat_PLU(), with m rows
  and m-1 columns
 
 invert transposed L and transposed Uu where Uu = U without zero row
 
 return result R as tmod_mat_single
 
 see comment on tmod_mat_solver_3arg for exact data layout
 '''
 cdef tmod_mat_t R
 cdef Py_ssize_t m=LU.matT[0].r
 tmod_mat_init_fast( R, m, m-1 )
 tmod_mat_solver_3arg( R, <long*>PD.array, LU.matT )
 return wrap_tmod_mat( R )

def export_Wti( tmod_mat_single WL_ti ):
 '''
 extract U'=Uu.T.I from tmod_mat_solver() result
 
 return U' as sage Matrix_integer_dense
 '''
 cdef Py_ssize_t i,j,i_plus, m=WL_ti.matT[0].r
 cdef Py_ssize_t cc=m-1
 cdef Matrix_integer_dense Uap=Matrix( cc, cc )
 cdef mp_limb_t* on_s_row
 cdef mpz_t*     on_t_row
 for i in range( cc ):
  i_plus=i+1
  on_s_row = WL_ti.matT[0].rows[i_plus]
  on_t_row = Uap._entries + i * cc
  for j in range(i_plus):
   mpz_set_ui( <mpz_ptr>on_t_row[j], on_s_row[j] )
 return Uap

def export_Lti( tmod_mat_single WL_ti ):
 '''
 extract L'=L.T.I from tmod_mat_solver() result
 
 WL_ti is virgin
 
 return L' as sage Matrix_integer_dense
 '''
 cdef Py_ssize_t i,j, m=WL_ti.matT[0].r
 cdef Py_ssize_t cc=m-1
 cdef Matrix_integer_dense Lap=Matrix( m, m )
 cdef Py_ssize_t delta_t=m+1
 cdef mp_limb_t* on_s_row=WL_ti.matT[0].entries
 cdef mpz_t*     on_t_row=Lap._entries
 for i in range( cc ):
  mpz_set_ui( <mpz_ptr>on_t_row[0], 1 )
  on_t_row += 1
  for j in range(cc-i):
   mpz_set_ui( <mpz_ptr>on_t_row[j], on_s_row[j] )
  on_s_row += m    # shift one row and one entry
  on_t_row += m    # shift one row
 mpz_set_ui( <mpz_ptr>on_t_row[0], 1 ) # last row has only 1
 return Lap

def clean_dict( d, s ):
 j=len(d)-s
 if j>0:
  kk,to_del=d.iterkeys(),[]
  for k in kk:
   to_del.append( k )
   j -= 1
   if not j:
    break 
  for k in to_del:
   del d[k]

def agnostic_array_export_big( agnostic_array A, long x, long y, mp_limb_t z ):
 '''
 treat A.array as array of mp_limb_t
 skip 1st x entries
 put other y entries into symmetric range, then take absolute value, call it u
 ignore all u < z
 return all bigger u as Python dict: u->count of occurencies, but make
  sure it contains not more than 100 keys (drop smaller u)
 if result is empty, return None
 '''
 cdef mp_limb_t* a=<mp_limb_t*>A.array
 a += x
 cdef long i
 cdef mp_limb_t u
 r=dict()
 for i in range(y):
  u=a[i]
  if u & 0x8000000000000000:
   u = -u
  if u<z:
   continue
  try:
   r[u] += 1
  except:
   r[u] = 1
   if len(r) > 100+10:
    clean_dict(r, 100 )
 if len(r) > 100:
  clean_dict(r, 100 )
 if len(r):
  return r
