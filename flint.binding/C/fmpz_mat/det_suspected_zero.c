// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>
#include "../mpz_square_mat/mpz_square_mat_.h"
#include "../nmod_mat/nmod_mat_.h"
#include "../fmpz/fmpz_.h"
#include "../ulong_extras/ulong_extras_.h"
#include "../fmpq/fmpq_.h"

/*
 gmp functions ..._ui() should work fine under Windoz/MPIR
 
 mpfr_get_uj() and mpfr_set_uj() might work, too 
*/

#define NDEBUG 0
/*
set NDEBUG to 0, to enable internal check in Dixon lifing algorithm
*/
#define LOUD_nmod_mat_in_det_divisor 1
#define LOUD_det_divisor_count_y 1

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

#define SQUARE(x,y,w) fmpz_get_mpz(y,w); mpz_mul(x,y,y)
int __inline__ static
log2_L2_fmpz_3arg(mpfr_t tgt,const fmpz* vec,slong n)
// 2*log2(L2 norm) rounded up. Return 1 on zero row
 {
  slong size=0,i,j;
  const fmpz* u;
  for(i=n,u=vec;i--;)
   {
    j=fmpz_size(u++);
    if(j>size)
     size=j;
   }
  if(0==size)
   return 1;
  j=size*FLINT_BITS;
  size = size*2*FLINT_BITS;
  mpz_t a,b,c; mpz_init2(a,size+FLINT_BITS); mpz_init2(b,j); mpz_init2(c,size);
  SQUARE( a, b, vec );
  for(i=n-1,u=vec+1;i--;u++)
   {
    SQUARE( c, b, u );
    mpz_add( a, a, c );
   }
  mpz_clear(c);
  mpfr_prec_t p=mpz_size(a)*FLINT_BITS;
  mpfr_t f; mpfr_init2(f,p);
  mpfr_set_z(f, a, MPFR_RNDU);
  mpfr_log2(tgt, f, MPFR_RNDU);
  mpfr_clear(f);
  mpz_clear(b); mpz_clear(a);
  return 0;
 }
#undef SQUARE

slong static __inline__
hadamard_3arg(mpfr_t b,const fmpz_mat_t m,mpfr_prec_t pr)
/*
upper bound on log2( 2*abs(m det) )
returns -1 if zero row found, smallest row index otherwise
*/
 {
  const slong n=m->r;
  slong smallest=0,j;
  // gcc warning: initialization from incompatible pointer type --- don't know
  //  how to fix
  const fmpz** const rows=m->rows;
  mpfr_t v; mpfr_init2(v,pr);
  mpfr_t u; mpfr_init2(u,pr);
  mpfr_t scratch; mpfr_init2(scratch,pr);
  #define HADAMARD_CLEAR mpfr_clear(scratch); mpfr_clear(u); mpfr_clear(v)
  if(log2_L2_fmpz_3arg( v, rows[0], n ))
   {
    HADAMARD_CLEAR;
    return -1;
   }
  mpfr_set(b, v, MPFR_RNDU);
  for(j=1;j<n;j++)
   {
    if(log2_L2_fmpz_3arg( u, rows[j], n ))
     {
      HADAMARD_CLEAR;
      return -1;
     }
    mpfr_add(b, b, u, MPFR_RNDU);
    if( mpfr_cmp(v, u)<0 )
     {
      smallest=j;
      mpfr_set(v, u, MPFR_RNDU);
     }
   }
  HADAMARD_CLEAR;
  mpfr_div_ui( b, b, 2, MPFR_RNDU ); // instead of taking root
  mpfr_add_ui( b, b, 1, MPFR_RNDU ); // instead of multiplying by 2
  #undef HADAMARD_CLEAR
  return smallest;
 }

#define SQUARE(x,y) mpz_mul(x,y,y)
void __inline__ static 
log2_L2_norm_4arg(mpfr_t tgt, mpz_square_mat_t A, slong k, slong n)
// log2(L2 norm) rounded down
 {
  mpz_ptr u;
  slong size=0,i,j;
  for(i=n,u=A->rows[k];i--;)
   {
    j=mpz_size(u++);
    if(j>size)
     size=j;
   }
  size = size*2*FLINT_BITS;
  mpz_t a,c; mpz_init2(a,size+FLINT_BITS); mpz_init2(c,size);
  u=A->rows[k];
  SQUARE( a, u );
  for(i=n-1,u++;i--;u++)
   {
    SQUARE( c, u );
    mpz_add( a, a, c );
   }
  mpz_clear(c);
  mpfr_prec_t p=mpz_size(a)*FLINT_BITS;
  mpfr_t f; mpfr_init2(f,p);
  mpfr_set_z(f, a, MPFR_RNDZ);
  mpfr_log2(tgt, f, MPFR_RNDZ);
  mpfr_div_ui( tgt, tgt, 2, MPFR_RNDZ );
  mpfr_clear(f);
  mpz_clear(a); 
 }
#undef SQUARE

mp_limb_t __inline__ static
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

slong __inline__ static 
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

void __inline__ static
det_divisor_init_b(mpz_ptr b,slong n,mpz_square_mat_t m)
 {
  slong i;
  mpz_ptr t;
  slong* s=m->mark;
  for(i=0,t=b;i<n;i++,t++)
   {
    mpz_init2(t, s[i]);
    mpz_set_si(t,  2*(i&1)-1);
   }
 }

void __inline__ static
det_divisor_clear_b(mpz_ptr b,slong n)
 {
  slong i;
  mpz_ptr t;
  for(i=n,t=b;i--;t++)
   mpz_clear(t);
  flint_free(b);
 }

mp_limb_t __inline__ static
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

void __inline__ static 
det_divisor_inverse_A(nmod_mat_t r, nmod_mat_t m, const fmpz_mat_t a, slong n)
 {
  memcpy( &r->mod, &m->mod, sizeof(m->mod) );
  nmod_mat_t s,t;
  nmod_mat_init_3arg(s, n, n);  memcpy( &s->mod, &m->mod, sizeof(m->mod) );
  nmod_mat_init_3arg(t, n, n);  memcpy( &t->mod, &m->mod, sizeof(m->mod) );
  fmpz_mat_get_nmod_mat(s,a);
  nmod_mat_inv(t,s);
  #if LOUD_nmod_mat_in_det_divisor
   flint_printf("A mod prime:\n");
   nmod_mat_print_pretty(s);
   flint_printf("\n\ninverted:\n");
   nmod_mat_print_pretty(t);
  #endif
  nmod_mat_transpose_square_tgt_virgin(r,t);
  #if LOUD_nmod_mat_in_det_divisor
   flint_printf("\n\ntransposed:\n");
   nmod_mat_print_pretty(r);
  #endif  
  nmod_mat_clear(t);
  nmod_mat_clear(s);
 }

void __inline__ static
print_limb_vector(char* m,const mp_limb_t* v,slong n)
 {
  flint_printf(m);
  slong i;
  for(i=0;i<n;i++)
   flint_printf("%wu ",v[i]);
  flint_printf("\n");
 }

void __inline__ static
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

void __inline__ static
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
    #if NDEBUG==0
     assert( 0 == mpz_tdiv_ui(pB, prime_p) );
    #endif
    mpz_divexact_ui(pB, pB, prime_p);
   }
 }

void __inline__ static
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

void __inline__ static
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
  det_divisor_clear_b(xA,n);
 }

void __inline__ static
det_divisor_ratnl_rcnstrction(mpz_t d,const nmod_mat_t y, slong k,
  slong n, mp_limb_t p, mp_limb_t log2_N, mp_limb_t log2_D
  #if NDEBUG==0
   ,mpz_square_mat_t A
  #endif
  )
 {
  slong i,Msize=FLINT_BITS*k;
  mpz_ptr zp,x_modulo_M=flint_malloc( sizeof(__mpz_struct)*n );
  //gmp_printf("x = ");
  for(i=0,zp=x_modulo_M;i<n;zp++,i++)
   {
    mpz_init2( zp, Msize );
    det_divisor_y_to_x( zp, i, y, k, n, p );
    //gmp_printf("%ZX ",zp);
   }
  //gmp_printf("\n");
  mpz_t M; mpz_init(M);
  mpz_ui_pow_ui(M,p,(mp_limb_t)k);                  // M=p**k
  #if NDEBUG==0
   //check that x_modulo_M*A equals original b modulo M
   det_divisor_xAbM_check(x_modulo_M,A,M,n);
  #endif
  det_divisor_rational_reconstruction(d, x_modulo_M, M, p, n, log2_N, log2_D);
  mpz_clear(M);
  det_divisor_clear_b( x_modulo_M, n );
 }

void __inline__ static 
fmpz_mat_det_divisor_8arg(mpz_t r,const fmpz_mat_t Ao, nmod_mat_t Amod,
  mpfr_t denominator_b, mpfr_prec_t pr, slong smallest_row, p_k_pk_t pp,
  mpz_t w)
// denominator_b >= log2(2*abs(A det))
 {
  mpz_square_mat_t A; mpz_square_mat_lazy_init_set(A, Ao);
  mp_limb_t denominator_b_i, numerator_b_i, bound;
  numerator_b_i=cramer_rule(denominator_b, A, pr, smallest_row);
  // numerator_b_i >= log2(numerator_b)
  denominator_b_i=mpfr_get_uj(denominator_b,MPFR_RNDU);
  if(denominator_b_i<FLINT_BITS+1)
   denominator_b_i=FLINT_BITS+1;
  bound=numerator_b_i+denominator_b_i;
  --denominator_b_i;
  // denominator_b_i >= log2(denominator_b)
  // bound >= log2(2*numerator_b*denominator_b) 
  slong max_i=dixon_lifting_max_i(bound, pp.p),i;
  // this many iterations Dixon lifting will take, i=0..max_i-1
  if(pp.k>1)
   {
    pp.k=1;
    init__p_pk__and__nmod(&pp, &Amod->mod);
   }
  const slong n=A->r;
  nmod_mat_t y_storage; nmod_mat_init_3arg(y_storage,max_i,n);
  nmod_mat_t Ainv;
  nmod_mat_init_3arg(Ainv,n,n); 
  det_divisor_inverse_A(Ainv, Amod, Ao, n);
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
    det_divisor_reduce_b(b_mod_p, b, n, pp.p);
    det_divisor_count_y(y_storage->rows[i], b_mod_p, Ainv, n);
    det_divisor_mul_add_divide(b, y_storage->rows[i], A_t, n, pp.p);
   }
  flint_free(b_mod_p);
  det_divisor_clear_b(b, n);
  mpz_square_mat_clear(A_t);
  nmod_mat_clear(Ainv);
  // for sanity check A is used
  det_divisor_ratnl_rcnstrction(r, y_storage, max_i, n, pp.p,
   numerator_b_i,denominator_b_i
   #if NDEBUG==0
    ,A
   #endif
  );
  gmp_printf("ratnl_rcnstrction() gave det divisor %ZX\n",r);
  nmod_mat_clear(y_storage);
  mpz_square_mat_clear(A);
  mpz_lcm(r, r, w);
 }

void __inline__ static
fmpz_mat_det_9arg(
  mpz_t tgt_det,
  const fmpz_mat_t A, nmod_mat_t Amod,
  mpfr_t hb, mpfr_prec_t pr, slong smallest_row,
  const p_k_pk_t pp, n_primes_rev_t it,mpz_t w)
/*
this subroutine only works when source matrice det is non-zero and a multiple
 of 2**126

w: known divisor of A determinant, w>1 
*/
 {
  mpz_t det_divisor; mpz_init(det_divisor);
  fmpz_mat_det_divisor_8arg( det_divisor, A, Amod, hb, pr, smallest_row, pp, w);
  gmp_printf("found det divisor %ZX\n",det_divisor);
  flint_printf("End of rails in fmpz_mat_det_8arg()\n");
  exit(1);
 }

mpfr_prec_t __inline__ static
mpfr_bitlength(mpfr_t s)
 {
  mpz_t B; mpz_init(B); mpfr_get_z(B, s, MPFR_RNDN);
  mpfr_prec_t q=mpz_sizeinbase(B,2);
  mpz_clear(B);
  return q;
 }

void __inline__ static
init_log2p_3arg(mpfr_t L,mpfr_t P,mpfr_prec_t q)
// initialize P then L
 {
  mpfr_init2(P, FLINT_BITS);
  mpfr_init2(L, q);
 }

void __inline__ static
det_suspected_zero_log2_w(mpfr_t r, const mpz_t w)
 {
  mpfr_prec_t p=mpz_size(w)*FLINT_BITS;
  mpfr_t t; mpfr_init2(t,p);
  mpfr_set_z(t,w,MPFR_RNDZ);
  mpfr_log2(r,t,MPFR_RNDZ);
  mpfr_clear(t);
 }

void __inline__ static
fmpz_mat_det_update_w(mpz_t u,n_primes_rev_t i,const mpz_t w)
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
  flint_rand_t rand_st; flint_randinit(rand_st);
  // add some randomness
  rand_st->__randval ^= (mp_limb_t)A->entries;
  // count Hadamard bound on det and remember which row is smallest
  mpfr_prec_t h_bits=hadamard_bits(A,rand_st);
  mpfr_t h_bound; mpfr_init2(h_bound,h_bits);
  slong smallest_row=hadamard_3arg(h_bound,A,h_bits);
  if(-1==smallest_row)
   {
    mpfr_clear(h_bound);
    mpz_set_ui(r,0);
    return;
   }
  // try to find prime p such that det A is non-zero modulo p
  //mpfr_printf("h_bits=%d h_bound=%10Rf\n",h_bits,h_bound);
  h_bits=mpfr_bitlength(h_bound)+5;
  mpfr_t prime_product; mpfr_init2(prime_product,h_bits);
  det_suspected_zero_log2_w(prime_product, W);
  const slong n=A->r;
  n_primes_rev_t it;
  mp_limb_t* scratch=NULL;
  p_k_pk_t pp; pp.p=0;
  nmod_mat_t Amod;
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
      pp.p=n_primes_rev_init(it, 0);
      scratch=flint_malloc( 4*(n-4)*sizeof(mp_limb_t) );
      init_log2p_3arg(log2_p, p, h_bits);
     }
    init__p_k_pk__and__nmod(&pp, &Amod->mod);
    fmpz_mat_get_nmod_mat(Amod, A); // TODO: speed-up this slow subroutine
    //flint_printf("selected p=%wu\n",pp.p);
    if( nmod_mat_det_mod_pk_4block(Amod,pp,scratch) )
     {
      mpz_t d; mpz_init_set(d, W);
      fmpz_mat_det_update_w(d, it, W);
      fmpz_mat_det_9arg(r, A, Amod, h_bound, h_bits, smallest_row, pp, it, d);
      mpz_clear(d);
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
  flint_randclear(rand_st);
 }

#undef NDEBUG
