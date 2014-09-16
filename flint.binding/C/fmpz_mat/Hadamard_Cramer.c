#include "fmpz_mat_.h"
#include "../mpz_square_mat/mpz_square_mat_.h"

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

mp_limb_t
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

slong
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
