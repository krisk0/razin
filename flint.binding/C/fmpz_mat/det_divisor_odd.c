// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>

#define po(x) flint_printf("%s\n",x);
#define SHOW_DIXON_RESULT 0
#define DIXON_INTERNAL_CHECK 0

void rational_reconstruction_2deg(mpz_t d,mpz_ptr x,slong n,mpz_t M,
  slong log2_M,mp_limb_t log2_N,mp_limb_t log2_D);
mp_limb_t tmod_mat_invert_transpose(tmod_mat_t R, const tmod_mat_t S);

#define STABLE_RP 0

#define THROW        \
 for(;;)                \
  {                        \
   m=gmp_urandomm_ui(rst,n); \
   if(0==b[m])                \
    break;                     \
  }                            \
 b[m]=2*(sign&1)-1;            \
 sign >>= 1;                   \
 ++thrown;

#if 0==STABLE_RP
 #include <time.h>
#endif

static __inline__ mp_limb_t
_20140914_form_b(mp_limb_t* b,slong n)
// throw 1 to 4 +1 or -1 into vector b (without collision)
 {
  // no good way to randomize in FLINT, using gmp procedures
  gmp_randstate_t rst; gmp_randinit_default(rst);
  #if 0==STABLE_RP
   gmp_randseed_ui(rst, (unsigned long)clock());
  #endif
  mp_limb_t sign=gmp_urandomb_ui(rst,4),thrown=0;
  slong i=13,k,m;
  THROW
  for(k=3;k--;)
   {
    if(  gmp_urandomm_ui(rst,23) < i  ) // 13/23 17/23 21/23
     break;
    THROW
    i += 4;
   }
  gmp_randclear(rst);
  #undef THROW
  #if SHOW_DIXON_RESULT
   flint_printf("0 b=");
   for(i=0;i<n;i++)
    flint_printf("%d,",(slong)b[i]);
   flint_printf("\n");
  #endif
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
  //mpfr_printf("log2 of %Rf equals %Rf\n",normF,tgt);
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
  //mpfr_printf("row0 L2 norm=%Rf\n",e);
  __mpfr_struct* eM=e;
  __mpfr_struct* eP;
  for(i=1,eP=e+1;i<n;i++,eP++)
   {
    pr=_20140914_log2_L2(eP, A->rows[i], n);
    //mpfr_printf("row %d L2 norm=%Rf\n",i,eP);
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
  MPFR_INIT2(h,b); mpfr_init2(c,b);
  slong i;
  __mpfr_struct* t;
  for(i=0,t=e; i<n; i++,t++)
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

static __inline__ void
_20140914_rp_Cramer(mp_limb_t* b, mpfr_t c,
  const mpfr_t h, slong n)
 {
  mp_limb_t b_norm=_20140914_form_b(b,n);
  {
   // add log2( b_norm )
   mpfr_t q,r;
   mpfr_init_set_ui(q,b_norm,MPFR_RNDU);
   mpfr_init(r);
   mpfr_log2(r,q,MPFR_RNDU);
   mpfr_add(c,c,r,MPFR_RNDU);
   mpfr_clear(r);
   mpfr_clear(q);
  }
  mpfr_div_ui(c,c,2,MPFR_RNDU);
  _20140914_lift_bound(c);
 }

slong _20140914_max_i(mpfr_t Ha, mpfr_t Cr)
 {
  #if DIXON_INTERNAL_CHECK
   mpfr_printf("hb=%Rf cb=%Rf\n",Ha,Cr);
  #endif
  return (mpfr_get_uj(Ha,MPFR_RNDU)+mpfr_get_uj(Cr,MPFR_RNDU)+
           FLINT_BITS-1)/FLINT_BITS;
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
   tr_neg->rows=flint_malloc(n*sizeof(void*));
   tr_neg->entries=t0;
   tr_neg->r=n;
   fmpz** const s_rows=s->rows;
   fmpz* g=s_rows[0];
   mp_limb_t* t1=x->entries;
   slong i,j;
   // transfer row 0, fill tr_neg->rows[]
   for(i=0; i<n; i++,g++,t0 += n,t1++)
    {
     tr_neg->rows[i]=t0;
     mpz_init(t0); fmpz_get_mpz(t0, g);
     t1[0]=mpz_to_t(t0);
    }
   for(j=1;j<n;j++)
    {
     // transfer row j
     g=s_rows[j];
     t0=tr_neg->entries+j;
     for(i=n;i--;g++,t0 += n,t1++)
      {
       mpz_init(t0); fmpz_get_mpz(t0, g);
       t1[0]=mpz_to_t(t0);
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
     sum += mpz_to_t(bP) * *q++;
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
    #if DIXON_INTERNAL_CHECK
     assert( 0==mpz_getlimbn(bP,0) );
    #endif
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
    #if DIXON_INTERNAL_CHECK
     assert( 0==mpz_getlimbn(bP,0) );
    #endif
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

static void
_20140914_check_x(mpz_ptr x,fmpz_mat_t a,mp_limb_t* b,mpz_t m,slong n)
// abort unless x*a-b is zero modulo m
 {
  mpz_ptr sI,s=flint_malloc(sizeof(__mpz_struct)*n);
  slong i,j;
  for(i=0,sI=s; i<n; i++,sI++)
   {
    mpz_init_set_si(sI, -(signed long) (b[i]) );
    for(j=0; j<n; j++)
     muladd_mpz_fmpz(sI, x+j, fmpz_mat_entry(a,j,i));
    mpz_mod(sI,sI,m);
   }
  j=0;
  for(i=0; i<n; i++)
   if(mpz_cmp_ui(s+i,0))
    {
     j=1;
     break;
    }
  if(0==j)
   return;
  flint_printf("x*A-B = ");
  for(i=0; i<n; i++)
   {
    if(i<n-1)
     gmp_printf("%ZX,",s+i);
    else
     gmp_printf("%ZX\n",s+i);
   }
  abort();
 }

static __inline__ void
_20140914_x_to_d(mpz_t d, mpz_ptr x, slong n, 
   slong bits,
   mp_limb_t hb, mp_limb_t cb
    #if DIXON_INTERNAL_CHECK
     ,fmpz_mat_t a, mp_limb_t* b
    #endif
  )
// feed correct args to rational_reconstruction_2deg()
 {
  mpz_t M; mpz_init_set_ui(M,1); mpz_mul_2exp(M,M,bits);
  #if SHOW_DIXON_RESULT
   gmp_printf("x=");
   int i;
   for(i=0; i<n; i++)
    gmp_printf("%Zd,",x+i);
   gmp_printf("\n");
  #endif
  #if DIXON_INTERNAL_CHECK
   gmp_printf("hb=%d cb=%d M=2**%d\n",hb,cb,bits);
   _20140914_check_x(x,a,b,M,n);
  #endif
  rational_reconstruction_2deg(d,x,n,M,bits,cb,hb);
  mpz_clear(M);
 }

static __inline__ void
_20140914_ratnl_rcnstrction(mpz_t r, const tmod_mat_t y, slong max_i, slong n,
   const mpfr_t hb, const mpfr_t cb
   #if DIXON_INTERNAL_CHECK
    ,fmpz_mat_t a, mp_limb_t* b
   #endif
  )
 {
  mpz_ptr xP,x=flint_malloc( sizeof(__mpz_struct)*n );
  slong i,p=max_i*FLINT_BITS;
  mp_limb_t hb_i=mpfr_get_uj(hb,MPFR_RNDU);
  mp_limb_t cb_i=mpfr_get_uj(cb,MPFR_RNDU);
  slong k=hb_i+cb_i;
  if(k==max_i*FLINT_BITS)
   k=-k;
  for(i=0,xP=x; i<n; xP++,i++)
   {
    mpz_init2(xP, p);
    _20140914_y_to_x(xP, i, y, max_i, n);
    if(k>0)
     mpz_mod_2x(xP, k);
   }
  _20140914_x_to_d(r, x, n, abs(k), hb_i-1, cb_i
    #if DIXON_INTERNAL_CHECK
     ,a, b
    #endif
  );
  clear_mpz_array(x, n);
 }

static void
_20140914_check_y0(mp_limb_t* y, mp_limb_t* b, fmpz_mat_t a, slong n)
// abort unless y*a equals b modulo T
 {
  slong i,j;
  mp_limb_t sum;
  for(i=0; i<n; i++)
   {
    sum=0;
    for(j=0; j<n; j++)
     sum += y[j] * fmpz_to_t( fmpz_mat_entry(a,j,i) );
    if( sum != b[i] )
     {
      gmp_printf("%Mu != %Mu\n",sum,b[i]);
      flint_printf("y*a not equals b, i=%d, y=",i);
      for(j=0; j<n; j++)
       gmp_printf("%Mu,",y[j]);
      flint_printf("\n");
      abort();
     }
   }
 }

static void
_20140914_check_yI(mp_limb_t* y, mpz_ptr b, fmpz_mat_t a, slong n)
// abort unless y*a equals b modulo T
 {
  slong i,j;
  mp_limb_t sum;
  for(i=0; i<n; i++)
   {
    sum=0;
    for(j=0; j<n; j++)
     sum += y[j] * fmpz_to_t( fmpz_mat_entry(a,j,i) );
    if( sum != mpz_to_t(b+i) )
     {
      flint_printf("y*a not equals b, i=%d\n",i);
      abort();
     }
   }
 }

static void
_20140914_check_A_inv_tr(tmod_mat_t inv_tr, fmpz_mat_t s, slong n)
// abort unless inv_tr * s transposed equals identity
 {
  mp_limb_t e=0;
  slong i,j;
  for(i=0; i<n; i++)
   for(j=0; j<n; j++)
    {
     mp_limb_t p=dot_modulo_t_kind0( inv_tr->rows[i], s->rows[j], n );
     if( (i==j) && (p!=1) )
      e=2;
     if( (i!=j) && (p!=0) )
      e=2;
     if(e)
      {
       gmp_printf("i/j=%d/%d p=0x%MX\n",i,j,p);
       gmp_printf("A_inv_tr check negative\n");
       abort();
      }
    }
 }


static __inline__ mp_limb_t
_20140914_Hadamard_bound(mpfr_t h, mpfr_t c, const fmpz_mat_t A, slong n)
/*
on entry h,c uninitialized

if 2*H.B. is larger than 2**(2*FLINT_BITS), 
 calculate h:=H.B.
 pre-calculate c:=C.B.
 return 0
else
 deinitialize h,c
 return log2( 2*H.B. ) rounded up
*/
 {
  __mpfr_struct* e=flint_malloc( sizeof(__mpfr_struct)*n );
  mpfr_prec_t b;
  slong smallest_row=_20140914_Hadamard(e,&b,A,n);
  MPFR_INIT2(h,b); mpfr_init2(c,b);
  slong i;
  __mpfr_struct* t;
  for(i=0,t=e; i<n; i++,t++)
   if(i!=smallest_row)
    mpfr_add(h,h,t,MPFR_RNDU);
  mpfr_set(c,h,MPFR_RNDU);   // c=sum log2 norm(i) excluding smallest_row
  mpfr_add(h,h,e+smallest_row,MPFR_RNDU);
  clear_mpfr_array(e,n);
  mpfr_div_ui(h,h,2,MPFR_RNDU);
  mpfr_add_ui(h,h,1,MPFR_RNDU);
  mp_limb_t hb_i=mpfr_get_uj(h,MPFR_RNDU);
  if(hb_i < 2*FLINT_BITS)
   {
    mpfr_clear(h);
    mpfr_clear(c);
    return hb_i;
   }
  return 0;
 }

mp_limb_t
fmpz_mat_det_divisor_odd(mpz_t r, mp_limb_t* det_mod_T, const fmpz_mat_t Ao)
/*
on entry r equals 1
 
returns log2( 2*H.B.(A) / r )  where r is the discovered divisor
*/
 {
  #if DIXON_INTERNAL_CHECK
   flint_printf("det_divisor_odd() source matrice\n"); fmpz_mat_print_pretty(Ao);
  #endif
  slong n=Ao->r;
  mp_limb_t hb_i;
  //rp_Hadamard_Cramer(Bo,hb,cb,Ao,n);
  mpfr_t hb,cb;
  hb_i=_20140914_Hadamard_bound(hb,cb,Ao,n);
  if(hb_i)
   {
    //  if H.B. is small, skip Dixon algorithm and leave r=1
    *det_mod_T=det_mod_t(Ao);
    return hb_i;
   }
  mp_limb_t* Bo=flint_calloc(sizeof(mp_limb_t),n);
  _20140914_rp_Cramer(Bo,cb,hb,n);
  // group 0: Bo, hb, cb
  slong max_i=_20140914_max_i(hb,cb);
  #if DIXON_INTERNAL_CHECK
   flint_printf("max_i=%d\n",max_i);
  #endif
  tmod_mat_t A_inv_tr;
  mpz_square_mat_t A_tr_neg;
  *det_mod_T=_20140914_matrix(A_inv_tr,A_tr_neg,Ao,n);
  #if DIXON_INTERNAL_CHECK
   tmod_mat_print_hex("A inv transposed modulo T:",A_inv_tr);
   _20140914_check_A_inv_tr(A_inv_tr,Ao,n);
  #endif
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
  #if DIXON_INTERNAL_CHECK
   _20140914_check_y0(Y->entries, Bo, Ao, n);
   po("y0 counted")
  #endif
  _20140914_update_B_step0( B, A_tr_neg, Y->entries, Bo, n );
  slong i;
  for(i=1;i<max_i;i++)
   {
    _20140914_count_y_main( Y->rows[i], B, A_inv_tr, n );
    #if DIXON_INTERNAL_CHECK
     _20140914_check_yI(Y->rows[i], B, Ao, n);
     po("yi counted")
    #endif
    _20140914_update_B_main( B, A_tr_neg, Y->rows[i], n );
   }
  // group 2
  clear_mpz_array(B, n);
  _20140914_ratnl_rcnstrction(r, Y, max_i, n, hb, cb
   #if DIXON_INTERNAL_CHECK
    ,Ao,Bo
   #endif
  );
  tmod_mat_clear(Y);
  // group 1
  mpz_square_mat_clear(A_tr_neg);
  tmod_mat_clear(A_inv_tr);
  // group 0
  mpfr_clear(cb);
  decrease_bound_mpz_2arg(hb,r);
  hb_i=mpfr_get_uj(hb,MPFR_RNDU);
  mpfr_clear(hb);
  flint_free(Bo);
  #if SHOW_DIXON_RESULT
   gmp_printf("fmpz_mat_det_divisor_odd() found divisor %Zd\n",r);
   gmp_printf("det modulo t=%Md\n",det_mod_T[0]);
   gmp_printf("fmpz_mat_det_divisor_odd() returning %Md\n",hb_i);
  #endif
  return hb_i;
 }

#undef STABLE_RP
#undef SHOW_DIXON_RESULT
#undef po
#undef DIXON_INTERNAL_CHECK
