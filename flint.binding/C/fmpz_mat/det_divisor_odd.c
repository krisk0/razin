// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>

void rational_reconstruction_2deg(mpz_t d,mpz_ptr x,slong n,mpz_t M,
  slong log2_M,mp_limb_t log2_N,mp_limb_t log2_D);
mp_limb_t tmod_mat_invert_transpose(tmod_mat_t R, const tmod_mat_t S);

#define STABLE_RP 1

#define THROW        \
 for(;;)                \
  {                       \
   m=n_randlimb(rst) % n;  \
   if(0==b[m])              \
    break;                   \
  }                           \
 b[m]=2*(sign&1)-1; \
 sign >>= 1;                   \
 ++thrown;

#if 0==STABLE_RP
 #include <time.h>
#endif

static __inline__ mp_limb_t
_20140914_form_b(mp_limb_t* b,slong n)
// throw 1 to 4 +1 or -1 into vector b (without collision)
 {
  flint_rand_t rst; flint_randinit(rst);
  // no good way to randomize in FLINT
  #if 0==STABLE_RP
   rst->__randval ^= (mp_limb_t)clock();
  #endif
  mp_limb_t sign=n_randlimb(rst),thrown=0;
  slong i=13,k,m;
  THROW
  for(k=3;k--;)
   {
    if(  n_randlimb(rst) % 23 < i  ) // 13/23 17/23 21/23
     break;
    THROW
    i += 4;
   }
  flint_randclear(rst);
  #undef THROW
  return thrown;
 }

__inline__ static mpfr_prec_t
_20140914_log2_L2(mpfr_t tgt,const fmpz* vec,slong n)
// 2*log2(L2 norm) rounded up
 {
  fmpz_t norm; fmpz_init(norm);
  slong i,i_is_big=0;
  square_L2_fmpz(norm,vec,n);
  i=fmpz_size(norm);
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
  if(i_is_big)
   return 1+FLINT_BITS;
  return mpfr_get_default_prec();
 }

static __inline__ slong
_20140914_Hadamard(__mpfr_struct* e,mpfr_prec_t* b,const fmpz_mat_t A,slong n)
 {
  slong smallest_row=0,i;
  mpfr_prec_t pr_max,pr;
  pr_max=pr=_20140914_log2_L2(e, A->rows[0], n);
  __mpfr_struct* eM=e;
  __mpfr_struct* eP;
  for(i=1,eP=e+1;i<n;i++,eP++)
   {
    pr=_20140914_log2_L2(eP, A->rows[i], n);
    if(pr>pr_max)
     pr_max=pr;
    if( mpfr_cmp(eP,eM) < 0 )
     {
      smallest_row=i;
      eM=eP;
     }
   }
  *b=pr_max;
  return smallest_row;
 }

static __inline__ void
_20140914_lift_bound(mpfr_t x)
 {
  if(mpfr_cmp_ui(x,FLINT_BITS)<0)
   mpfr_set_ui(x,FLINT_BITS,MPFR_RNDU);
 }

static __inline__ void
Hadamard_Cramer(mpfr_t h,mpfr_t c,mp_limb_t rp_norm,const fmpz_mat_t A,slong n)
 {
  __mpfr_struct* e=flint_malloc( sizeof(__mpfr_struct)*n );
  mpfr_prec_t b;
  slong smallest_row=_20140914_Hadamard(e,&b,A,n);
  mpfr_init2(h,b); mpfr_init2(c,b);
  slong i;
  __mpfr_struct* t;
  for(i=0,t=e;i<n;i++,t++)
   if(i!=smallest_row)
    mpfr_add(h,h,t,MPFR_RNDU);
  {
   mpfr_t q,r;
   mpfr_init_set_ui(q,rp_norm,MPFR_RNDU);
   mpfr_init(r);
   mpfr_log2(r,q,MPFR_RNDU);
   mpfr_set(c,r,MPFR_RNDU);
   mpfr_clear(r);
   mpfr_clear(q);
  }
  mpfr_add(c,c,h,MPFR_RNDU);
  mpfr_div_ui(c,c,2,MPFR_RNDU);
  mpfr_add(h,h,e+smallest_row,MPFR_RNDU);
  mpfr_div_ui(h,h,2,MPFR_RNDU);
  _20140914_lift_bound(h);
  mpfr_add_ui(h,h,1,MPFR_RNDU);
  clear_mpfr_array(e,n);
 }

static __inline__ void
rp_Hadamard_Cramer(mp_limb_t*b, mpfr_t h, mpfr_t c, fmpz_mat_t A, slong n)
/*
h>=log2(2*H.B.(A)), c>=log2(Cramer bound)
h>=1+FLINT_BITS, c>=FLINT_BITS
*/
 {
  mp_limb_t b_norm=_20140914_form_b(b,n);
  Hadamard_Cramer(h,c,b_norm,A,n);
  _20140914_lift_bound(c);
 }

slong _20140914_max_i(mpfr_t Ha, mpfr_t Cr)
 {
  mpfr_t t; mpfr_init2(t, mpfr_get_prec(Ha));
  mpfr_add(t,Ha,Cr,MPFR_RNDU);
  mpfr_div_ui(t,t,FLINT_BITS,MPFR_RNDU);
  slong r=mpfr_get_si(t,MPFR_RNDU);
  mpfr_clear(t);
  return r;
 }

static __inline__ mp_limb_t
_20140914_matrix(tmod_mat_t inv_tr,mpz_square_mat_t tr_neg,const fmpz_mat_t s,
  slong n)
 {
  tmod_mat_t x;
  tmod_mat_init_fast(inv_tr,n,n);
  tmod_mat_init_fast(x,n,n);
  // transpose  s -> tr_neg  and  let x = s mod 2**64
  // mpz_square_mat_init_transpose_fmpz(tr_neg,s);
  {
   mpz_ptr t0=flint_malloc(n*n*sizeof(__mpz_struct));
   tr_neg->entries=t0;
   tr_neg->r=n;
   fmpz** const s_rows=s->rows;
   fmpz* g=s_rows[0];
   mp_limb_t* t1=x->entries;
   slong i,j;
   // transfer row 0, fill tr_neg->rows[]
   for(i=0;i<n;i++,g++,t0 += n,t1++)
    {
     tr_neg->rows[i]=t0;
     mpz_init(t0); fmpz_get_mpz(t0, g);
     t1[0]=flint_mpz_get_ui(t0);
    }
   for(j=1;j<n;j++)
    {
     // transfer row j
     g=s_rows[j];
     t0=tr_neg->entries+j;
     for(i=n;i--;g++,t0 += n,t1++)
      {
       mpz_init(t0); fmpz_get_mpz(t0, g);
       t1[0]=flint_mpz_get_ui(t0);
      }
    }
  }
  // tr_neg->mark uninitialized, except ->mark[]
  mp_limb_t det=tmod_mat_invert_transpose(inv_tr,x);
  tmod_mat_clear(x);
  mpz_square_mat_negate(tr_neg);
  mpz_square_mat_mark_biggest(tr_neg);
  return det;
 }

static __inline__ void
_20140914_count_y_step0(mp_limb_t* y,mp_limb_t* b,tmod_mat_t a,slong n)
/*
 y = b*a transposed

 a virgin
*/ 
 {
  mp_limb_t* q=a->entries;
  slong i,j;
  for(i=0;i<n;i++)
   {
    mp_limb_t sum=0;
    for(j=0;j<n;j++)
     sum += b[j] * *q++;
    y[i]=sum;
   }
 }

static __inline__ void
_20140914_count_y_main(mp_limb_t* y,mpz_ptr b,tmod_mat_t a,slong n)
 {
  mp_limb_t* q=a->entries;
  slong i,j;
  for(i=0;i<n;i++)
   {
    mp_limb_t sum=0;
    mpz_ptr bP=b;
    for(j=n;j--;bP++)
     sum += mpz_mod_T(bP) * *q++;
    y[i]=sum;
   }
 }

static __inline__ void
_20140914_update_B_step0(mpz_ptr b, const mpz_square_mat_t a, mp_limb_t* y, 
  mp_limb_t* b0, slong n)
/*
b=(b0+y*a transposed) / 2**64

a virgin
*/
 {
  mpz_ptr bP=b;
  mpz_ptr aP=a->entries;
  slong i,j;
  for(i=0;i<n;i++,bP++)
   {
    mpz_set_si(bP, (unsigned long) b0[i]);
    for(j=0;j<n;j++,aP++)
     flint_mpz_addmul_ui(bP,aP,y[j]);
    assert( 0==mpz_getlimbn(bP,0) );
    mpz_shift_right_1limb(bP);
   }
 }

static __inline__ void
_20140914_update_B_main(mpz_ptr b, const mpz_square_mat_t a, mp_limb_t* y, 
  slong n)
/*
b=(b+y*a transposed) / 2**64

a virgin
*/
 {
  mpz_ptr bP=b;
  mpz_ptr aP=a->entries;
  slong i,j;
  for(i=0;i<n;i++,bP++)
   {
    for(j=0;j<n;j++,aP++)
     flint_mpz_addmul_ui(bP,aP,y[j]);
    assert( 0==mpz_getlimbn(bP,0) );
    mpz_shift_right_1limb(bP);
   }
 }

static __inline__ void
_20140914_y_to_x(mpz_t x,slong i,const tmod_mat_t y, slong y_rows, slong y_cols)
 {
  slong j=y_rows-1;
  mp_limb_t const * p=y->rows[j]+i;
  while( (0==p[0]) && (j>=0) )
   {
    --j;
    p -= y_cols;
   }
  if(j<0)
   {
    x->_mp_size=0;
    return;
   }
  x->_mp_size=j+1;
  mp_limb_t* q=x->_mp_d+j;
  for(;j--;p -= y_cols)
   *q-- = *p;
  *q = *p;
 }

static __inline__ void
_20140914_x_to_d(mpz_t d, mpz_ptr x, slong n, slong k, 
  const mpfr_t hb_plus_1, const mpfr_t cb)
// feed correct args to rational_reconstruction_2deg()
 {
  mpz_t M; mpz_init_set_ui(M,1); mpz_mul_2exp(M,M,k *= FLINT_BITS);
  rational_reconstruction_2deg(d,x,n,M,k,
    mpfr_get_uj(hb_plus_1,MPFR_RNDU)-1,
    mpfr_get_uj(cb,MPFR_RNDU)
   );
  mpz_clear(M);
 }

static __inline__ void
_20140914_ratnl_rcnstrction(mpz_t r, const tmod_mat_t y, slong max_i, slong n,
  const mpfr_t hb,const mpfr_t cb)
 {
  mpz_ptr xP,x=flint_malloc( sizeof(__mpz_struct)*n );
  slong i;
  for(i=0,xP=x;i<n;xP++,i++)
   {
    mpz_init2(xP, max_i);
    _20140914_y_to_x(xP, i, y, max_i, n);
   }
  _20140914_x_to_d(r, x, n, max_i, hb, cb);
  clear_mpz_array(x, n);
 }

mp_limb_t
fmpz_mat_det_divisor_odd(mpz_t r, mp_limb_t* det_mod_T, fmpz_mat_t Ao)
/*
returns log2( 2*H.B.(A) / r )  where r is the discovered divisor
*/
 {
  slong n=Ao->r;
  mp_limb_t* Bo=flint_calloc(sizeof(mp_limb_t),n);
  mpfr_t hb,cb;
  rp_Hadamard_Cramer(Bo,hb,cb,Ao,n);
  // group 0: Bo, hb, cb
  slong max_i=_20140914_max_i(hb,cb);
  tmod_mat_t A_inv_tr;
  mpz_square_mat_t A_tr_neg;
  *det_mod_T=_20140914_matrix(A_inv_tr,A_tr_neg,Ao,n);
  // group 1: A_inv_tr, A_tr_neg
  tmod_mat_t Y; tmod_mat_init_fast(Y, max_i, n);
  mpz_ptr B=flint_malloc( sizeof(__mpz_struct)*n );
  det_divisor_init_zero_b(B, n, A_tr_neg);
  // group 2: Y, B
  #if 0
   for i in range(max_i):
    Yi = B*A inverted modulo T
    B = (B+Yi*(-A)) / T  --- must divide exactly
   x = sum Yi*T**i
   rational reconstruction(x)
  #endif
  _20140914_count_y_step0( Y->entries, Bo, A_inv_tr, n );
  _20140914_update_B_step0( B, A_tr_neg, Y->entries, Bo, n );
  slong i;
  for(i=1;i<max_i;i++)
   {
    _20140914_count_y_main( Y->rows[i], B, A_inv_tr, n );
    _20140914_update_B_main( B, A_tr_neg, Y->rows[i], n );
   }
  // group 2
  clear_mpz_array(B, n);
  _20140914_ratnl_rcnstrction(r, Y, max_i, n, hb, cb);
  tmod_mat_clear(Y);
  // group 1
  mpz_square_mat_clear(A_tr_neg);
  tmod_mat_clear(A_inv_tr);
  // group 0
  mpfr_clear(cb);
  decrease_bound_mpz_2arg(hb,r);
  mp_limb_t hb_i=mpfr_get_uj(hb,MPFR_RNDU);
  mpfr_clear(hb);
  flint_free(Bo);
  return hb_i;
 }

#undef STABLE_RP
