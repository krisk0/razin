// This program is part of RAZIN
// Copyright Денис Крыськов 2014

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include <flint/fmpz_mat.h>
#include "../fmpz/fmpz_.h"
#include "../ulong_extras/ulong_extras_.h"

void fmpz_mat_det_5arg(mpz_t r,const fmpz_mat_t A,mpfr_t h_bound,
 slong smallest_row,mp_limb_t good_p);

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
  i=(slong)mpfr_get_ui(r, D);             // convert result to int
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

slong static
_hadamard(mpfr_t b,const fmpz_mat_t m,mpfr_prec_t pr)
// upper bound on log2( m det )
 {
  #define D MPFR_RNDU
  #define L2( tgt, p )       \
   {                             \
    SQUARE( tgt, p[0] )             \
    for(i=1;i<n;i++)                 \
     {                                \
      SQUARE( scratch, p[i] );         \
      mpfr_add( tgt, tgt, scratch, D ); \
     }                                  \
    mpfr_log2( tgt, tgt, D );          \
   }
  #define SQUARE( s_tgt, f )           \
   fmpz_get_mpfr_macro(s_tgt, f, D);  \
   mpfr_mul(s_tgt, s_tgt, s_tgt, D);
  const slong n=m->r;
  slong smallest=0,i,j;
  const fmpz** const rows=m->rows;
  mpfr_t v; mpfr_init2(v,pr);
  mpfr_t u; mpfr_init2(u,pr);
  mpfr_t scratch; mpfr_init2(scratch,pr);
  L2( v, rows[0] );
  mpfr_set(b, v, D);
  for(j=1;j<n;j++)
   {
    L2( u, rows[j] );
    mpfr_add(b, b, u, D);
    if( mpfr_cmp(v,u)<0 )
     {
      smallest=j;
      mpfr_set(v,u,D);
     }
   }
  mpfr_div_ui( b, b, 2, D );
  mpfr_add_ui( b, b, 1, D );
  return smallest;
  #undef SQUARE
  #undef L2
  #undef D
 }

void
fmpz_mat_det_suspected_zero(mpz_t r,const fmpz_mat_t A,const mpz_t W)
/*
 det A is known to be a multiple of W and suspected to be zero, W>2**64,
  A dimension > 4

 count det A
*/
 {
  flint_rand_t rand_st; flint_randinit(rand_st);
  // add some randomness
  rand_st->__randval ^= (mp_limb_t)A->entries;
  // count Hadamard bound on det and remember which row is smallest
  mpfr_prec_t h_bits=hadamard_bits(A,rand_st);
  mpfr_t h_bound; mpfr_init2(h_bound,h_bits);
  slong smallest_row=_hadamard(h_bound,A,h_bits);
  // find prime p such that det A is non-zero modulo p, or make sure that it is
  //  zero
  //mpfr_printf("h_bits=%d h_bound=%10Rf\n",h_bits,h_bound);
  mpfr_t prime_product; mpfr_init2(prime_product,h_bits);
  mpfr_set_z(prime_product, W, MPFR_RNDZ);
  mpfr_log2(prime_product, prime_product, MPFR_RNDZ);
  const slong n=A->r;
  n_primes_rev_t it;
  mp_limb_t* scratch;
  p_k_pk_t pp; pp.p=0;
  nmod_mat_t Amod;
  mpfr_t log2_p;
  while( mpfr_cmp(prime_product,h_bound) < 0 )
   {
    if(pp.p)
     {
      pp.p=n_primes_rev_next(it);
      //pp.p==1 when range exhausted, not checking it
     }
    else
     {
      nmod_mat_init_square_2arg(Amod,n);
      pp.p=n_primes_rev_init(it,0);
      scratch=flint_malloc( 4*(n-4)*sizeof(mp_limb_t) );
      mpfr_init2(log2_p, h_bits);
     }
    init__p_k_pk__and__nmod(&pp,&Amod->mod);
    fmpz_mat_get_nmod_mat(Amod, A);
    //flint_printf("selected p=%wu\n",pp.p);
    if( nmod_mat_det_mod_pk_4block(Amod,pp,scratch) )
     {
      fmpz_mat_det_5arg(r, A, h_bound, smallest_row, pp.p);
      pp.p=UWORD_MAX;
      break;
     }
    mpfr_set_ui(log2_p, pp.p_deg_k, MPFR_RNDZ);
    mpfr_log2(log2_p, log2_p, MPFR_RNDZ);
    mpfr_add(prime_product, prime_product, log2_p, MPFR_RNDZ);
   }
  if(pp.p)
   {
    mpfr_clear(log2_p);
    flint_free(scratch);
    n_primes_rev_clear(it);
    nmod_mat_clear(Amod);
   }
  if(pp.p != UWORD_MAX)
   mpz_set_ui(r,0);
  mpfr_clear(prime_product);
  mpfr_clear(h_bound);
 }
