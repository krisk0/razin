// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>
#include "../mpz_square_mat/mpz_square_mat_.h"
#include "../fmpz/fmpz_.h"
#include "../ulong_extras/ulong_extras_.h"
#include "../ulong_extras/profile_.h"
#include "../fmpq/fmpq_.h"
#include "fmpz_mat_.h"
#include "../nmod_mat/nmod_mat_.h"

/*
 gmp functions ..._ui() should work fine under Windoz/MPIR
 
 mpfr_get_uj() and mpfr_set_uj() might work, too 
*/

// no asserts in code below, if checks disabled
#define LOUD_nmod_mat_in_det_divisor 0
#define LOUD_det_divisor_count_y 0
#define DIXON_INTERNAL_CHECK 0
#define LOUD_DET_RESULT 0

mpfr_prec_t static
hadamard_bits(const fmpz_mat_t m,flint_rand_t r_st)
// sample some entries of m and guess bit-size of Hadamard bound on det: 
//  log2(det) = 1/2*n*log2(n)*log2( squared entry ) =
//                  n*log2(n)*log2(         entry ) 
 {
  #define RAND_ROW n_randlimb(r_st) % n
  #define D MPFR_RNDN
  #define SAMPLE(z)                                       \
   fmpz_get_mpfr_macro(z, m->rows[RAND_ROW][RAND_ROW], D); \
   (void)mpfr_abs(z, z, D);
  slong i,n=m->r;
  mpfr_t r; mpfr_init(r);
  mpfr_t nn; mpfr_init(nn);
  SAMPLE( r );
  for(i=9;i--;)
   {
    SAMPLE( nn )
    if(mpfr_cmp_ui(nn,0))                    // skip zero
     {
      (void)mpfr_log2(nn, nn, D );
      (void)mpfr_add(r, r, nn, D);
     }
   }
  (void)mpfr_div_ui(r, r, 10, D);            // r=average log2(m entry)
  (void)mpfr_set_ui( nn, (unsigned long)m->r, D );
  (void)mpfr_mul(r, r, nn, D );              // multiply by n
  mpfr_log2( nn, nn, D );
  (void)mpfr_mul(r, r, nn, D);              // multiply by log2(n)
  i=(slong)mpfr_get_ui(r, D)+10;          // convert result to int, add 10
  mpfr_clear(nn); 
  mpfr_clear(r);
  n=(slong)mpfr_get_default_prec();   // 53 on 64-bit Linux
  if(i<n)
   i=n;                            // not less than default precision
  return (mpfr_prec_t)i;
  #undef SAMPLE
  #undef D
  #undef RAND_ROW
 }

__inline__ static void
square_L2_fmpz(fmpz_t r,const fmpz* vec,slong n)
// on entry r=0. on exit r=squared L2 norm of the vector 
 {
  fmpz_mul(r,vec,vec);
  const fmpz* u;
  fmpz_t x; fmpz_init(x);
  for(u=vec+1,--n;n--;u++)
   {
    fmpz_mul(x,u,u);
    fmpz_add(r,r,x);
   }
  fmpz_clear(x);
 }

__inline__ static void
square_L2_mpz(mpz_t r,mpz_srcptr vec,slong n)
// on entry r=0. on exit r=squared L2 norm of the vector 
 {
  mpz_mul(r,vec,vec);
  mpz_srcptr u;
  mpz_t x; mpz_init(x);
  for(u=vec+1,--n;n--;u++)
   {
    mpz_mul(x,u,u);
    mpz_add(r,r,x);
   }
  mpz_clear(x);
 }

__inline__ static int
log2_L2_fmpz_3arg(mpfr_t tgt,const fmpz* vec,slong n)
// 2*log2(L2 norm) rounded up
// if row is zero, return 1 and do not initialize tgt
// else return 0 and initialize tgt to a suitable precision
 {
  fmpz_t norm; fmpz_init(norm);
  slong i,i_is_big=0;
  square_L2_fmpz(norm,vec,n);
  i=fmpz_size(norm);
  if( 0==i )
   {
    fmpz_clear(norm);
    return 1;
   }
  if(i>2)
   {
    i=2;
    i_is_big=1;
   }
  mpfr_t normF; mpfr_init2(normF,i*FLINT_BITS);
  fmpz_get_mpfr(normF,norm,MPFR_RNDU);
  fmpz_clear(norm);
  if(i_is_big)
   mpfr_init2(tgt,1+FLINT_BITS);
  else
   mpfr_init(tgt);
  mpfr_log2(tgt, normF, MPFR_RNDU);
  mpfr_clear(normF);
  return 0;
 }

static __inline__ void
mpfr_copy_bound(mpfr_t r,const mpfr_t s)
 {
  mpfr_init2(r, mpfr_get_prec(s));
  mpfr_set(r, s, MPFR_RNDU);
 }

static __inline__ void
mpfr_add_bound(mpfr_t r,const mpfr_t s)
 {
  mpfr_prec_t s_prec=mpfr_get_prec(s);
  if(s_prec > mpfr_get_prec(r) )
   {
    mpfr_t n; mpfr_init2(n,s_prec);
    mpfr_swap(n,r);
    mpfr_add(r,n,s,MPFR_RNDU);
    mpfr_clear(n);
   }
  else
   mpfr_add(r,r,s,MPFR_RNDU);
 }

static __inline__ slong
hadamard_2arg(mpfr_t b,const fmpz_mat_t m)
/*
upper bound on log2( 2*abs(m det) )
returns -1 if zero row found, smallest row index otherwise

b on entry is uninitialized
b on exit is initialized iff no zero row found
*/
 {
  const slong n=m->r;
  slong smallest=0,j;
  // gcc warning: initialization from incompatible pointer type --- don't know
  //  how to fix
  const fmpz** const rows=m->rows;
  mpfr_t v,u;
  if(log2_L2_fmpz_3arg( v, rows[0], n ))
   return -1;
  mpfr_copy_bound(b, v);
  // v and b must be freed
  for(j=1;j<n;j++)
   {
    if(log2_L2_fmpz_3arg( u, rows[j], n ))
     {
      mpfr_clear(b); mpfr_clear(v);
      return -1;
     }
    mpfr_add_bound(b, u);
    if( mpfr_cmp(u, v)<0 )
     {
      smallest=j;
      mpfr_swap(v, u);
     }
    mpfr_clear(u);
    // v and b must be freed
   }
  mpfr_clear(v);
  mpfr_div_ui( b, b, 2, MPFR_RNDU ); // instead of taking root
  mpfr_add_ui( b, b, 1, MPFR_RNDU ); // instead of multiplying by 2
  return smallest;
 }

__inline__ static void
log2_L2_norm_4arg(mpfr_t tgt, const mpz_square_mat_t A, slong k, slong n)
// log2(L2 norm) rounded down. tgt initialized by caller
 {
  mpz_t norm; mpz_init(norm);
  slong i;
  square_L2_mpz(norm,A->rows[k],n);
  i=mpz_size(norm);
  if(i>2)
   i=2;
  mpfr_t normF; mpfr_init2(normF,i*FLINT_BITS);
  mpfr_set_z(normF,norm,MPFR_RNDU);
  mpz_clear(norm);
  mpfr_log2(tgt, normF, MPFR_RNDU);
  mpfr_div_ui(tgt,tgt,2,MPFR_RNDU);
  mpfr_clear(normF);
 }

__inline__ static mp_limb_t
cramer_rule(const mpfr_t den_bound, 
  mpz_square_mat_t A, mpfr_prec_t pr, slong k)
// returns log2(Cramer bound on numerator) rounded up
 {
  mp_limb_t bI;
  const slong n=A->r;
  mpfr_t u,v,w;
  mpfr_init2(w,pr);
  mpfr_init2(u,pr);
  log2_L2_norm_4arg(u, A, k, n);
  //mpfr_printf("k=%d row norm=%10Rf\n",k,u);
  mpfr_sub(w, den_bound, u, MPFR_RNDU);  // w=den_bound / min row norm
  mpfr_sub_ui(u, w, 1, MPFR_RNDU);       // u=den_bound/2/min row norm
  // vector norm = square root of n
  mpfr_init(v);
  mpfr_set_ui(v, n, MPFR_RNDU);
  mpfr_log2(w, v, MPFR_RNDU);            // w=b norm * 2
  mpfr_clear(v);
  mpfr_div_ui(w, w, 2, MPFR_RNDU);       // w=b norm
  mpfr_add(w, w, u, MPFR_RNDU);          // w=b norm*(...)
  mpfr_clear(u);
  bI=mpfr_get_uj(w, MPFR_RNDU);
  mpfr_clear(w);
  if(bI<FLINT_BITS)
   bI=FLINT_BITS;
  return bI;
 }

__inline__ static slong
dixon_lifting_max_i(mp_limb_t b,mp_limb_t p)
 {
  slong r;
  mpfr_t i; mpfr_init(i);
  mpfr_t pF; mpfr_init2(pF,FLINT_BITS);
  mpfr_t bF; mpfr_init2(bF,FLINT_BITS);
  // p might be big, using mpfr_set_uj()
  (void)mpfr_set_uj(pF,p,MPFR_RNDZ);
  (void)mpfr_log2(pF, pF, MPFR_RNDZ);  // pF <= log2(p)
  (void)mpfr_set_uj(bF,b,MPFR_RNDU);   // bF >= b
  (void)mpfr_div(i, bF, pF, MPFR_RNDU); // i >= b/log2(p)
  r=(slong)mpfr_get_uj(i,MPFR_RNDU);   // r >= i
  mpfr_clear(bF);
  mpfr_clear(pF);
  mpfr_clear(i);
  return r;
 }

__inline__ static mp_limb_t
det_divisor_reduce_b(mp_limb_t* tgt, mpz_srcptr src, slong n,mp_limb_t p)
 {
  slong i;
  mpz_srcptr s;
  for(i=0,s=src;i<n;i++,s++)
   {
    if( mpz_cmp_ui(s,0) >= 0)
     tgt[i] = mpz_tdiv_ui(s,p);   //  mpz_tdiv_ui returns absolute value 
    else
     tgt[i] = p-mpz_tdiv_ui(s,p); //   of remainder, so negate it
   }
  #if LOUD_det_divisor_count_y
   flint_printf("ZZ b=");
   for(i=0;i<n;i++)
    gmp_printf("%Zd ",src+i);
   flint_printf("\n");
  #endif
 }

__inline__ static void
det_divisor_inverse_A(nmod_mat_t r,const p_k_pk_t const* pp,nmod_mat_t m, 
  const fmpz_mat_t a, slong n)
 {
  #if LOUD_nmod_mat_in_det_divisor
   int check_result;
  #endif
  p_k_pk_t ppM;
  if(pp->k==1)
   memcpy( &r->mod, &m->mod, sizeof(m->mod) );
  else
   {
    // k changed, must re-calculate mod 
    ppM.p=ppM.p_deg_k=pp->p;
    ppM.k=1;
    init_nmod_from_pp(&r->mod, &ppM);
   }
  // r->mod calculated, now for the entries
  nmod_mat_t t; nmod_mat_init_3arg(t, n, n);

  // must modify m->mod, then restore it
  nmod_t s; memcpy( &s, &m->mod, sizeof(s) );
  if(pp->k>1)
   // we want inverse modulo p, not p**k
   nmod_init( &m->mod, pp->p );
  else
   mod_flint_style( &m->mod, pp );
  memcpy( &t->mod, &m->mod, sizeof(s) );
  nmod_mat_inv(t,m);
  memcpy( &m->mod, &s, sizeof(s) );
  
  // check inverse
  #if LOUD_nmod_mat_in_det_divisor
   nmod_mat_mul_pk_classical(r, t, m);
   check_result=nmod_mat_is_one(r);
  #endif
  // check inverse positive iff check_result=1
  
  #if LOUD_nmod_mat_in_det_divisor
   flint_printf("A mod prime:\n");
   nmod_mat_print_pretty(s);
   flint_printf("\n\ninverted:\n");
   nmod_mat_print_pretty(t);
   if( !check_result )
    {
     flint_printf("det_divisor_inverse_A(): bad inverse\n");
     abort();
    }
  #endif
  nmod_mat_transpose_square_tgt_virgin(r,t);
  #if LOUD_nmod_mat_in_det_divisor
   flint_printf("\n\ntransposed:\n");
   nmod_mat_print_pretty(r);
  #endif  
  nmod_mat_clear(t);
 }

__inline__ static void
print_limb_vector(char* m,const mp_limb_t* v,slong n)
 {
  flint_printf(m);
  slong i;
  for(i=0;i<n;i++)
   flint_printf("%wu ",v[i]);
  flint_printf("\n");
 }

__inline__ static void
det_divisor_count_y(mp_limb_t* y,const mp_limb_t* b,const nmod_mat_t m,slong dim)
// y := b*transposed m
 {
  const mp_limb_t** const m_row=m->rows;
  const mp_limb_t n=m->mod.n;
  const mp_limb_t nI=m->mod.ninv;
  const mp_limb_t norm=m->mod.norm;
  slong i,j;
  #if defined(VECTOR_DOT_TAIL) && defined(SPEEDUP_NMOD_RED3)
  for(i=0;i<dim;i++)
   {
    const mp_limb_t* c=m_row[i];
    VECTOR_DOT_HEAD( b[0], c[0] );
    for(j=1;j<dim;j++)
     VECTOR_DOT_BODY( b[j], c[j] );
    VECTOR_DOT_TAIL( y[i], n,nI,norm );
   }
  #else
   #error cant compile this
  #endif
  #if LOUD_det_divisor_count_y
   print_limb_vector("b=", b, dim);
   print_limb_vector("y=", y, dim);
  #endif
 }

__inline__ static void
det_divisor_mul_add_divide(mpz_ptr b,const mp_limb_t* y,
  const mpz_square_mat_t m,slong n,mp_limb_t prime_p)
// b := (b+y*transposed m)/prime_p, b[i] allocated to correct size already  
 {
  slong i,j;
  mpz_ptr pB,pM;
  const mp_limb_t* pY;
  slong* w=m->mark;
  const mpz_ptr* m_rows=m->rows;
  for(i=0,pB=b;i<n;i++,pB++)
   {
    mpz_t s; mpz_init2(s, w[i]);
    // TODO: what if mp_limb_t is not unsigned long?
    for(j=n,pM=m_rows[i],pY=y;j--;pM++,pY++)
     {
      mpz_mul_ui(s, pM, *pY);
      mpz_add(pB, pB, s);
     }
    mpz_clear(s);
    // b[i] should divide by p exactly
    #if DIXON_INTERNAL_CHECK
     assert( 0 == mpz_tdiv_ui(pB, prime_p) );
    #endif
    mpz_divexact_ui(pB, pB, prime_p);
   }
 }

__inline__ static void
det_divisor_y_to_x(mpz_t xI, slong i, const nmod_mat_t y, slong y_rows, slong y_cols,
  mp_limb_t q)
 {
  slong j=y_rows-1;
  mp_limb_t const * p=y->rows[j]+i;
  mpz_set_ui(xI,*p);
  // TODO: what if mp_limb_t is not unsigned long?
  for(;j--;)
   {
    p -= y_cols;
    mpz_mul_ui(xI,xI,q);
    mpz_add_ui(xI,xI,*p);
   }
 }

__inline__ static void
det_divisor_xAbM_check(mpz_ptr x,mpz_square_mat_t A,const mpz_t M,slong n)
 {
  mpz_ptr zp,xA=flint_malloc( sizeof(__mpz_struct)*n );
  slong i;
  for(i=0,zp=xA;i<n;zp++,i++)
   mpz_init(zp);
  mpz_square_mat_mul_vec_mat_modulo(xA,x,A,M);
  for(i=0,zp=xA;i<n;zp++,i++)
   {
    //mpz_add_si(zp,-((i&1) ? 1 : -1)); // this addition should give zero 
    gmp_printf("xA[%d]=%ZX\n",i,zp);
    if( i&1 )
     mpz_sub_ui(zp, zp, 1);
    else
     mpz_add_ui(zp, zp, 1);
    assert( (0==mpz_cmp_ui(zp,0) || (0==mpz_cmp(zp,M)) ) );
   }
  flint_printf("xA=b modulo M  check positive\n");
  clear_mpz_array(xA,n);
 }

__inline__ static void
det_divisor_ratnl_rcnstrction(mpz_t d,const nmod_mat_t y, slong k,
  slong n, mp_limb_t p, mp_limb_t log2_N, mp_limb_t log2_D
  #if DIXON_INTERNAL_CHECK
   ,mpz_square_mat_t A
  #endif
  )
 {
  slong i,Msize=FLINT_BITS*k;
  mpz_ptr zp,x_modulo_M=flint_malloc( sizeof(__mpz_struct)*n );
  for(i=0,zp=x_modulo_M;i<n;zp++,i++)
   {
    mpz_init2( zp, Msize );
    det_divisor_y_to_x( zp, i, y, k, n, p );
   }
  mpz_t M; mpz_init(M);
  mpz_ui_pow_ui(M,p,(mp_limb_t)k);                  // M=p**k
  #if DIXON_INTERNAL_CHECK
   gmp_printf("M = %ZX\n",M);
   //check that x_modulo_M*A equals original b modulo M
   det_divisor_xAbM_check(x_modulo_M,A,M,n);
  #endif
  det_divisor_rational_reconstruction(d, x_modulo_M, M, p, n, log2_N, log2_D);
  mpz_clear(M);
  clear_mpz_array( x_modulo_M, n );
 }

void
print_y_by_A(mp_limb_t* y, const nmod_mat_t A, slong n)
 {
  mp_limb_t* x=flint_malloc( n*sizeof(mp_limb_t) );
  nmod_mat_mul_vec_left(x, y, A);
  flint_printf("y*A mod p=");
  slong i;
  for(i=0;i<n;i++)
   flint_printf("%wu ",x[i]);
  flint_printf("\n");
  flint_free(x);
 }

void
test_y_by_A(mp_limb_t* y, const nmod_mat_t A, slong n,const mp_limb_t* b)
 {
  slong i,ok=1;
  mp_limb_t* x=flint_malloc( n*sizeof(mp_limb_t) );
  nmod_mat_mul_vec_left(x, y, A);
  for(i=n;i--;)
   if(x[i] != b[i])
    ok=0;
  if( !ok )
   {
    flint_printf("A=\n");   nmod_mat_print_pretty(A); flint_printf("\n\n");
    flint_printf("y*A mod p=");
    for(i=0;i<n;i++)
     flint_printf("%wu ",x[i]);
    flint_printf("\ny*A not equals b\n");
    flint_printf("A->mod: %wu, %wu, %wu\n",A->mod.n,A->mod.ninv,A->mod.norm);
    abort();
   }
  flint_free(x);
 }

__inline__ static void
fmpz_mat_det_divisor_7arg(mpz_t r,const fmpz_mat_t Ao, nmod_mat_t Amod,
  mpfr_t denominator_b, mpfr_prec_t pr, slong smallest_row, p_k_pk_t pp
  #if RAT_REC_TAKES_D_SERIOUSLY==0  
   ,mpz_t w
  #endif
  )
// denominator_b >= log2(2*abs(A det))
 {
  #if RAT_REC_TAKES_D_SERIOUSLY
   // divide denominator bound by known det divisor r
   decrease_bound_fmpz(denominator_b,pr,r);
   #if LOUD_DET_BOUND
    mpfr_printf("H.B. for Dixon (log2): %Rf\n",denominator_b);
   #endif
  #endif
  mpz_square_mat_t A; mpz_square_mat_lazy_init_set(A, Ao);
  mp_limb_t denominator_b_i, numerator_b_i, bound;
  numerator_b_i=cramer_rule(denominator_b, A, pr, smallest_row);
  // numerator_b_i >= log2(numerator_b)
  denominator_b_i=mpfr_get_uj(denominator_b,MPFR_RNDU);
  if(denominator_b_i<FLINT_BITS+1)
   denominator_b_i=FLINT_BITS+1;
  bound=numerator_b_i+denominator_b_i;
  #if LOUD_DET_BOUND
   flint_printf("log2 bound for Dixon: %llX\n",bound);
  #endif
  --denominator_b_i;
  // denominator_b_i >= log2(denominator_b)
  // bound >= log2(2*numerator_b*denominator_b) 
  slong max_i=dixon_lifting_max_i(bound, pp.p),i;
  // this many iterations Dixon lifting will take, i=0..max_i-1
  const slong n=A->r;
  nmod_mat_t y_storage; nmod_mat_init_3arg(y_storage,max_i,n);
  nmod_mat_t Ainv;
  nmod_mat_init_3arg(Ainv,n,n); 
  det_divisor_inverse_A(Ainv, &pp, Amod, Ao, n);
  mpz_square_mat_t A_t; mpz_square_mat_init_transpose(A_t,A);
  mpz_square_mat_mark_biggest(A_t);
  mpz_square_mat_negate(A_t);
  mpz_ptr b=flint_malloc( sizeof(__mpz_struct)*n );
  det_divisor_init_b(b, n, A_t);
  mp_limb_t* b_mod_p=(mp_limb_t*)flint_malloc(sizeof(mp_limb_t)*n);
  // TODO: pre-reduce b
  #if 0
   b={-1,1,...} --- integer vector
   for i in range(max_i):
    y := b*A inverted modulo p
    store y into row no. i of y storage
    b := (b-y*A)/p  --- must divide exactly
   x := sum yI*p**i
   rational reconstruction(x)
  #endif
  for(i=0;i<max_i;i++)
   {
    #if LOUD_det_divisor_count_y
     flint_printf("Dixon loop, i=%d\n",i);
    #endif
    det_divisor_reduce_b(b_mod_p, b, n, pp.p);
    det_divisor_count_y(y_storage->rows[i], b_mod_p, Ainv, n);
    #if LOUD_det_divisor_count_y
     //print_y_by_A(y_storage->rows[i], Amod, n);
     test_y_by_A(y_storage->rows[i], Amod, n, b_mod_p);
    #endif
    det_divisor_mul_add_divide(b, y_storage->rows[i], A_t, n, pp.p);
   }
  flint_free(b_mod_p);
  clear_mpz_array(b, n);
  mpz_square_mat_clear(A_t);
  nmod_mat_clear(Ainv);
  // for sanity check A is used
  det_divisor_ratnl_rcnstrction(r, y_storage, max_i, n, pp.p,
   numerator_b_i,denominator_b_i
   #if DIXON_INTERNAL_CHECK
    ,A
   #endif
  );
  //gmp_printf("ratnl_rcnstrction() gave det divisor %ZX\n",r);
  nmod_mat_clear(y_storage);
  mpz_square_mat_clear(A);
  #if RAT_REC_TAKES_D_SERIOUSLY==0  
   mpz_lcm(r, r, w);
  #endif
 }

__inline__ static void
fmpz_mat_det_9arg(
  mpz_t tgt_det,
  const fmpz_mat_t A, nmod_mat_t Amod,
  mpfr_t hb, mpfr_prec_t pr, slong smallest_row,
  p_k_pk_t pp, n_primes_rev_t it,mp_limb_t xmod
  #if RAT_REC_TAKES_D_SERIOUSLY==0  
   ,mpz_t w
  #endif
  )
/*
this subroutine only works when source matrice det is non-zero and a multiple
 of 2**126

w: known divisor of A determinant, w>1 
*/
 {
  #if LOUD_DET_RESULT
   gmp_printf("initial det divisor %ZX\n",
    #if RAT_REC_TAKES_D_SERIOUSLY
     tgt_det
    #else
     w
    #endif
    );
  #endif
  // Do not subtract log2(det divisor) twice --- keep original bound
  mpfr_t hb0; mpfr_copy_bound(hb0, hb);
  MARK_TIME(t0);
  fmpz_mat_det_divisor_7arg(tgt_det, A, Amod, hb, pr, smallest_row, pp
    #if RAT_REC_TAKES_D_SERIOUSLY==0  
     , w
    #endif
   );
  #if LOUD_DET_RESULT
   gmp_printf("after rat. rec.: %ZX\n",tgt_det);
  #endif
  DUMP_TIME("fmpz_mat_det_divisor_7arg()",t0);
  MARK_TIME(t1);
  fmpz_mat_det_modular_given_divisor_8arg(tgt_det,Amod,hb0,pr,&pp,it,xmod,A);
  #if LOUD_DET_RESULT
   gmp_printf("fmpz_mat_det_9arg() result: %ZX\n",tgt_det);
  #endif
  DUMP_TIME("fmpz_mat_det_modular_given_divisor_8arg()",t1);
  mpfr_clear(hb0);
 }

__inline__ static void
init_log2p_3arg(mpfr_t L,mpfr_t P,mpfr_prec_t q)
// initialize P then L
 {
  mpfr_init2(P, FLINT_BITS);
  mpfr_init2(L, q);
 }

void __inline__ static
det_suspected_zero_log2_w(mpfr_t r, const mpz_t w)
 {
  mpfr_t t; mpfr_init2(t,mpz_size(w)*FLINT_BITS);
  mpfr_set_z(t,w,MPFR_RNDZ);
  mpfr_log2(r,t,MPFR_RNDZ);
  mpfr_clear(t);
 }

void __inline__ static
fmpz_mat_det_update_w(mpz_t u,n_primes_rev_t i)
// re-iterate thru primes again, multiply w by those primes, skipping last
 {
  mp_limb_t last_p=n_primes_rev_show_again(i),curr_p;
  n_primes_rev_reset(i);
  // TODO: what if mp_limb_t is not unsigned long?
  while(1)
   {
    curr_p=n_primes_rev_show_again(i);
    if(curr_p == last_p)
     break;                  // i returned to where it was, go away
    mpz_mul_ui(u,u,curr_p);
    (void)n_primes_rev_next(i);
   }
 }

void
fmpz_mat_det_suspected_zero(mpz_t r,const fmpz_mat_t A,const mpz_t W)
/*
 det A is known to be a multiple of W and suspected to be zero,
  W>1, A dimension > 4

 count det A
*/
 {
  MARK_TIME(t_st)
  // count Hadamard bound on det and remember which row is smallest
  mpfr_t h_bound; // hadamard_2arg() initialises h_bound, 
                  // fmpz_mat_det_suspected_zero() releases it
  MARK_TIME(t_h)
  slong smallest_row=hadamard_2arg(h_bound,A);
  DUMP_TIME("hadamard_2arg()",t_h)
  if(-1==smallest_row)
   {
    mpz_set_ui(r,0);
    return;
   }
  // try to find prime p such that det A is non-zero modulo p
  //mpfr_printf("h_bits=%d h_bound=%10Rf\n",h_bits,h_bound);
  mpfr_prec_t h_bits=mpfr_get_prec(h_bound);
  mpfr_t prime_product; mpfr_init2(prime_product,h_bits);
  det_suspected_zero_log2_w(prime_product, W);
  const slong n=A->r;
  n_primes_rev_t it;
  mp_limb_t* scratch=NULL;
  p_k_pk_t pp; pp.p=0;
  nmod_mat_t Amod,Amod_ori;
  mpfr_t log2_p,p;
  while( mpfr_cmp(prime_product,h_bound) < 0 )
   {
    if(pp.p)
     {
      pp.p=n_primes_rev_next(it);
      //pp.p==1 when range exhausted, not checking it
     }
    else
     {
      nmod_mat_init_square_2arg(Amod, n);
      nmod_mat_init_square_2arg(Amod_ori, n);
      pp.p=n_primes_rev_init(it, 0);
      scratch=flint_malloc( 4*(n-4)*sizeof(mp_limb_t) );
      init_log2p_3arg(log2_p, p, h_bits);
     }
    init__p_k_pk__and__nmod(&pp, &Amod->mod);
    fmpz_mat_get_nmod_mat(Amod, A); // TODO: speed-up this slow subroutine
    //flint_printf("selected p=%wu\n",pp.p);
    nmod_mat_copy_entries(Amod_ori,Amod);
    mp_limb_t xmod=nmod_mat_det_mod_pk_4block(Amod,pp,scratch);
    if( xmod )
     {
      memcpy( &Amod_ori->mod, &Amod->mod, sizeof(Amod->mod) );
      #if RAT_REC_TAKES_D_SERIOUSLY
       mpz_set(r, W);
       fmpz_mat_det_update_w(r, it);
       fmpz_mat_det_9arg(r, A, Amod_ori, h_bound, h_bits, smallest_row, pp, it,
        xmod);
      #else
       mpz_t d; mpz_init_set(d, W);
       fmpz_mat_det_update_w(d, it);
       fmpz_mat_det_9arg(r, A, Amod_ori, h_bound, h_bits, smallest_row, pp, it,
        xmod, d);
       mpz_clear(d);
      #endif
      pp.p=UWORD_MAX;
      break;
     }
    mpfr_set_uj(p, pp.p_deg_k, MPFR_RNDZ);
    mpfr_log2(log2_p, p, MPFR_RNDZ);
    mpfr_add(prime_product, prime_product, log2_p, MPFR_RNDZ);
   }
  if(scratch)
   {
    mpfr_clear(log2_p);
    mpfr_clear(p);
    flint_free(scratch);
    n_primes_rev_clear(it);
    nmod_mat_clear(Amod);
   }
  if(pp.p != UWORD_MAX)
   mpz_set_ui(r,0);
  mpfr_clear(prime_product);
  mpfr_clear(h_bound);
  DUMP_TIME("count_det_suspected_zero()",t_st)
 }

#undef NDEBUG
