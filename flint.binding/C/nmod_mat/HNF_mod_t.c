#if WANT_ASSERT_IN_HNF_MOD_T
 #include <assert.h>
 #define ASSERT(x) assert(x)
#else
 #define ASSERT(x)
#endif

#define EMPTY_CACHE_CELL ( (mp_limb_t)-1 )
#define SENIOR_BIT_MASK ( UWORD(1) << (FLINT_BITS-1) )
#define HEAVY_CACHE_SIZE 4

#define LOUD_ROW_OP 0

void
_27022015_print_diag_and_above(nmod_mat_t d,slong clean_column,mp_limb_t m)
 {
  slong n=d->r,i,j;
  nmod_mat_t c; nmod_mat_init_3arg(c,n,n); memset(&c->mod,0,sizeof(c->mod));
  if(m)
   c->mod.n=m;
  for(i=n;i--;)
   for(j=n;j--;)
    {
     if( (j>i) && (i<=clean_column) )
      nmod_mat_entry(c,j,i) = 0;
     else
      {
       if(0==m)
        nmod_mat_entry(c,j,i) = nmod_mat_entry(d,j,i);
       else
        nmod_mat_entry(c,j,i) = nmod_mat_entry(d,j,i) % m;
      }
    }
  nmod_mat_print_pretty(c);
  nmod_mat_clear(c);
 }

static __inline__ void
_27022015_col_cache(nmod_mat_t A,mp_limb_t* cache,slong row,slong col)
 {
  mp_limb_t aIJ=nmod_mat_entry(A,row,col);
  if(0==aIJ)
   {
    cache[1]=EMPTY_CACHE_CELL;
    return;
   }
  if( aIJ & SENIOR_BIT_MASK )
   {
    cache[1]=-aIJ;
    cache[0]=(row<<1)^1;
   }
  else
   {
    cache[1]= aIJ;
    cache[0]= row<<1   ;
   }
 }

static __inline__ void
_27022015_bubble_sort(mp_limb_t* a,slong siz_minus_1)
 {
  while(siz_minus_1)
   {
    slong i;
    for(i=0;i<siz_minus_1;i++)
     if( a[2*i+1] > a[2*i+3] )
      {
       mp_limb_t t;
       t=a[2*i+1]; a[2*i+1]=a[2*i+3]; a[2*i+3]=t;
       t=a[2*i  ]; a[2*i  ]=a[2*i+2]; a[2*i+2]=t;
      }
    siz_minus_1--;
   }
 }

static __inline__ slong
_27022015_row_op_full(nmod_mat_t A,mp_limb_t* sc,
  slong x_row,mp_limb_t u,mp_limb_t w,
  slong y_row,mp_limb_t v,mp_limb_t z,
  mp_limb_t g,slong col,slong size)
/*
General case: 
 x*u + y*v replaces x
 x*w - y*z replaces y
*/
 {
  #if LOUD_ROW_OP
   gmp_printf("finally u, v / w, z = %Mu, %Mu / %Mu, %Mu\n",u,v,w,z);
  #endif
  mp_limb_t* xP=A->rows[x_row]+col;
  mp_limb_t* yP=A->rows[y_row]+col;
  --size;
  ASSERT( xP[0]*u + yP[0]*v == g);
  ASSERT( xP[0]*w - yP[0]*z == 0);
  *xP++ = g;
  *yP++ = 0;
  tmod_vec_mul(       sc, xP, size, u);
  tmod_vec_muladd(    sc, yP, size, v); // s=x*u+y*v
  tmod_vec_mul_3arg(  yP, size, -z );
  tmod_vec_muladd(    yP, xP, size, w );
  memcpy(             xP, sc, size*sizeof(mp_limb_t) );
  return x_row;
 }

static __inline__ slong
_27022015_row_op(nmod_mat_t A,mp_limb_t* scrtch,
  slong alpha_row, mp_limb_t alpha_num, slong alpha_mul,
  slong betta_row, mp_limb_t betta_num, slong betta_mul,
  slong col)
/*
-1 if cache exhausted
*/
 {
  if(betta_num==EMPTY_CACHE_CELL)
   return -1;
  mp_limb_t u,v,g;
  g=n_gcd_ui_positive_4arg(&u,&v,alpha_num,betta_num);
  #if LOUD_ROW_OP
   gmp_printf("alpha/row betta/row | g u v = %Mu/%d %Mu/%d | %Mu %Mu %Mu\n",
    alpha_num,alpha_row,betta_num,betta_row,g,u,v);
   gmp_printf("alpha/betta multiplyiers: %d %d\n",alpha_mul,betta_mul);
  #endif
  if(alpha_mul)
   u=-u;
  if(betta_mul)
   v=-v;
  // A[alpha_row,col]*u + A[betta_row,col]*v = g
  slong size=A->r-col;
  if(0==v)
   {
    if(alpha_mul)
     tmod_vec_neg_2arg( A->rows[alpha_row]+col, size );
    // betta_num is a multiple of alpha_num
    v=betta_num / alpha_num;
    if(!betta_mul)
     v=-v;
    ASSERT( nmod_mat_entry(A,alpha_row,col)*(-v) == 
            nmod_mat_entry(A,betta_row,col) );
    #if LOUD_ROW_OP
     gmp_printf("adding to row %d row %d multiplied by %Mu\n",
      betta_row,alpha_row,v);
    #endif
    tmod_vec_muladd( A->rows[betta_row]+col, A->rows[alpha_row]+col, v, size );
    return alpha_row;
   }
  if(0==u)
   {
    if(betta_mul)
     tmod_vec_neg_2arg( A->rows[betta_row]+col, size );
    // alpha_num is a multiple of betta_num
    u=alpha_num / betta_num;
    if(!alpha_mul)
     u=-u;
    ASSERT( nmod_mat_entry(A,betta_row,col)*(-u) == 
            nmod_mat_entry(A,alpha_row,col) );
    #if LOUD_ROW_OP
     gmp_printf("adding to row %d row %d multiplied by %Mu\n",
      alpha_row,betta_row,u);
    #endif
    tmod_vec_muladd( A->rows[alpha_row]+col, A->rows[betta_row]+col, u, size );
    return betta_row;
   }
  alpha_mul=1-2*alpha_mul;
  betta_mul=1-2*betta_mul;
  ASSERT( alpha_mul*alpha_num == nmod_mat_entry(A,alpha_row,col) );
  ASSERT( betta_mul*betta_num == nmod_mat_entry(A,betta_row,col) );
  if(alpha_row<betta_row)
   {
    return _27022015_row_op_full(A,scrtch,
     alpha_row,u,betta_num/g*alpha_mul,
     betta_row,v,alpha_num/g*betta_mul,
     g,col,size);
   }
  return _27022015_row_op_full(A,scrtch,
   betta_row,v,alpha_num/g*betta_mul,
   alpha_row,u,betta_num/g*alpha_mul,
   g,col,size);
 }

static __inline__ slong
_27022015_row_op_3arg(nmod_mat_t A,mp_limb_t* scrtch,mp_limb_t* cache,slong col)
// engage rows as indicated by 2 initial cache cells
 {
  return _27022015_row_op(A,scrtch,
   cache[0]>>1, cache[1], cache[0]&1,
   cache[2]>>1, cache[3], cache[2]&1,
   col);
 }

static __inline__ slong
_27022015_row_op_4arg(nmod_mat_t A,mp_limb_t* scrtch,mp_limb_t* cache,
  slong alpha_row,slong col)
// engage row alpha_row to row from cache
 {
  return _27022015_row_op(A,scrtch,
   alpha_row, nmod_mat_entry(A,alpha_row,col), 0,
   cache[0]>>1,                 cache[1],      cache[0]&1,
   col);
 }

static __inline__ slong
_27022015_fix_diag(nmod_mat_t A,slong i)
//make sure A[i,i] is a degree of 2
 {
  mp_limb_t aII=nmod_mat_entry(A,i,i);
  ASSERT( aII );
  if(!is_degree_of_2(aII))
   {
    mp_limb_t s;
    // zap trailing zeroes
    count_trailing_zeros(s,aII);
    s=aII>>s;
    // invert modulo t
    s=t_invmod(s);
    // multiply row by s
    mp_limb_t* p=&nmod_mat_entry(A,i,i);
    tmod_vec_mul_3arg(p,A->r-i,s);
    ASSERT( is_degree_of_2( p[0] ) );
   }
  return 0;
 }

static __inline__ int
_27022015_finalize_row_op(nmod_mat_t A,slong last_row,slong col)
/*
if necessary, switch rows
if necessary, reduce diagonal element so it is a degree of 2
return 0
*/
 {
  if(last_row != col)
   MP_PTR_SWAP( A->rows[col], A->rows[last_row] );
  return _27022015_fix_diag(A,col);
 }

static __inline__ int
_27022015_diag_and_down_light(nmod_mat_t A,mp_limb_t* scrtch,
  slong i,slong row_len)
 {
  mp_limb_t col_cache[6];
  _27022015_col_cache( A, col_cache+0,i,i );
  _27022015_col_cache( A, col_cache+2*1,i+1,i );
  if(row_len==2)
   _27022015_col_cache( A, col_cache+2*2,i+2,i );
  _27022015_bubble_sort( col_cache, row_len );
  if(col_cache[1]==EMPTY_CACHE_CELL)
   return 1;
  slong curr_row=_27022015_row_op_3arg(A,scrtch,col_cache,i);
  if(curr_row<0)
   {
    // cache exhausted, clean-up and go away
    return _27022015_finalize_row_op(A,col_cache[0]>>1,i);
   }
  slong prev_row=curr_row;
  if(row_len==2)
   {
    curr_row=_27022015_row_op_4arg(A,scrtch,col_cache+4,curr_row,i);
    if(curr_row>=0)
     prev_row=curr_row;
   }
  return _27022015_finalize_row_op(A,prev_row,i);
 }

static __inline__ slong
_27022015_find_odd(const nmod_mat_t A,slong col)
 {
  slong i,max_plus_1=A->r;
  for(i=col;i<max_plus_1;i++)
   {
    if( 1 & nmod_mat_entry(A,i,col) )
     return i;
   }
  return -1;
 }

static __inline__ void
_27022015_make_1(mp_limb_t* rowP,slong size_minus_1)
 {
  ASSERT( 1 & rowP[0] );
  mp_limb_t q=t_invmod(rowP[0]);
  *rowP++ = 1;
  tmod_vec_mul_3arg(rowP,size_minus_1,q);
 }

static __inline__ void
_27022015_zap_with_1(nmod_mat_t A,mp_limb_t* pivot_row,slong col,slong row_no,
   slong size_minus_1)
 {
  mp_limb_t* this_row=A->rows[row_no]+col;
  mp_limb_t aJI=this_row[0];
  if(aJI)
   tmod_vec_muladd(this_row+1,pivot_row,size_minus_1,-aJI);
 }

static __inline__ int
_27022015_find_odd_invert_zap(nmod_mat_t A,slong col,slong row_len)
/*
Try and find invertible element among A[col,col],A[col+1,col],...
use it as pivot
return 0 on failure
*/
 {
  slong i=_27022015_find_odd(A,col);
  if(i<0)
   return 0;
  if( i != col )
   MP_PTR_SWAP( A->rows[col], A->rows[i] );
  mp_limb_t* pivot_P=A->rows[col]+col;
  _27022015_make_1(pivot_P,row_len);
  ++pivot_P;
  for(i=col+1;i<A->r;i++)
   _27022015_zap_with_1(A,pivot_P,col,i,row_len);
  return 1;
 }

static __inline__ void
_27022015_maybe_cache(mp_limb_t* cache,slong cache_size,slong row,mp_limb_t e)
 {
  mp_limb_t c0=0;
  if(e&SENIOR_BIT_MASK)
   {
    c0=1;
    e = -e;
   }
  if(e>cache[2*cache_size-1])
   return;
  cache[2*cache_size-1]=e;
  cache[2*cache_size-2]=(row<<1)^c0;
  for(row=cache_size-1;row;row--)
   {
    #if 0
     gmp_printf("Comparing key at %d which is %MX agains key at %d which is"
      " %MX\n",row,cache[2*row+1],row-1,cache[2*row-1]);
    #endif
    if(cache[2*row+1] < cache[2*row-1])
     {
      c0=cache[2*row+1]; cache[2*row+1]=cache[2*row-1]; cache[2*row-1]=c0;
      c0=cache[2*row  ]; cache[2*row  ]=cache[2*row-2]; cache[2*row-2]=c0;
     }
    else
     break;
   }
 }

static __inline__ void
_27022015_heavy_cache(mp_limb_t* cache,slong cache_size,nmod_mat_t A,slong col)
 {
  slong i;
  for(i=col;i<A->r;i++)
   {
    mp_limb_t aIJ=nmod_mat_entry(A,i,col);
    if(aIJ)
     _27022015_maybe_cache(cache,cache_size,i,aIJ);
   }
 }

static __inline__ int
_27022015_diag_and_down_heavy(nmod_mat_t A,mp_limb_t* scrtch,slong col)
 {
  // select 4 best rows
  mp_limb_t cache[2*HEAVY_CACHE_SIZE];
  cache[1]=cache[3]=cache[5]=cache[7]=EMPTY_CACHE_CELL;
  _27022015_heavy_cache(cache,HEAVY_CACHE_SIZE,A,col);
  #if LOUD_ROW_OP
   gmp_printf("_diag_and_down_heavy(): col=%d  cache keys: %MX %MX %MX %MX\n",
    col,cache[1],cache[3],cache[5],cache[7]);
  #endif  
  if(cache[1]==EMPTY_CACHE_CELL)
   return 1;
  // engage 2 champions
  slong curr_row=_27022015_row_op_3arg(A,scrtch,cache,col);
  #if LOUD_ROW_OP
   gmp_printf("_diag_and_down_heavy(): col=%d champions fight result: %d\n",
    col,curr_row);
  #endif  
  if(curr_row<0)
   return _27022015_finalize_row_op(A,cache[0]>>1,col);
  slong last_row=curr_row,i;
  // 2 more encounters
  for(i=2;i<HEAVY_CACHE_SIZE;i++)
   {
    curr_row=_27022015_row_op_4arg(A,scrtch,cache+2*i,curr_row,col);
    #if LOUD_ROW_OP
     gmp_printf("diag_and_down_heavy() done with cache cell %d, curr_row=%d\n",
      i,curr_row);
    #endif
    if(curr_row>=0)
     last_row=curr_row;
    else
     break;
   }
  if(last_row != col)
   MP_PTR_SWAP( A->rows[col], A->rows[last_row] );
  if(curr_row<0)
   return _27022015_fix_diag(A,col);
  #if LOUD_ROW_OP
   gmp_printf("diag_and_down_heavy(): afterparty starts\n");
   _27022015_print_diag_and_above(A,col,0);
  #endif
  for(i=col+1;i<A->r;i++)
   {
    mp_limb_t aIJ=nmod_mat_entry(A,i,col);
    if(aIJ)
     {
      #if LOUD_ROW_OP
       gmp_printf("diag_and_down_heavy(): afterparty, i=%d\n",i);
      #endif
      slong mul=0;
      if( aIJ & SENIOR_BIT_MASK )
       {
        mul=1;
        aIJ=-aIJ;
       }
      mp_limb_t pivot=nmod_mat_entry(A,col,col);
      _27022015_row_op(A,scrtch, col,pivot,0, i,aIJ,mul, col);
      #if LOUD_ROW_OP
       gmp_printf("diag_and_down_heavy(): done with i=%d\n",i);
       _27022015_print_diag_and_above(A,col-1,0);
      #endif
     }
   }
  return _27022015_fix_diag(A,col);
 }

static __inline__ int
_27022015_diag_and_down(nmod_mat_t A,slong i,mp_limb_t* scrtch)
 {
  slong row_len=A->r-i-1;
  if(row_len)
   {
    if( _27022015_find_odd_invert_zap(A,i,row_len) )
     return 0;
    if(row_len>2)
     return _27022015_diag_and_down_heavy(A,scrtch,i);
    return _27022015_diag_and_down_light(A,scrtch,i,row_len);
   }
  else
   {
    mp_limb_t* aIIp=A->rows[i]+i;
    mp_limb_t aII=aIIp[0];
    if(!aII)
     return 1;
    #if LOUD_ROW_OP
     gmp_printf("lower-right before fixing: %Mu\n",aII);
    #endif
    if(!is_degree_of_2(aII))
     {
      count_trailing_zeros( row_len, aII );
      #if LOUD_ROW_OP
       gmp_printf("zeroes: %d\n",row_len);
      #endif
      aIIp[0]=UWORD(1)<<row_len;
     }
    return 0;
   }
 }

static __inline__ void
HNF_mod_t_det(fmpz_t r,const nmod_mat_t m)
 {
  fmpz_set_ui(r,UWORD(1));
  mp_limb_t** rows=m->rows;
  ulong k=0,j;
  for(j=m->r;j--;)
   {
    mp_limb_t mJJ=rows[j][j];
    if(mJJ)
     fmpz_mul_ui(r, r, mJJ);
    else
     ++k;
   }
  if(k)
   fmpz_mul_2exp(r, r, k*FLINT_BITS);
 }

static __inline__ void
HNF_mod_t_det_no_zeros(fmpz_t r,const nmod_mat_t m)
 {
  mp_limb_t** rows=m->rows;
  slong j=m->r-1;
  fmpz_set_ui(r,rows[j][j]);
  for(;j--;)
   fmpz_mul_ui(r, r, rows[j][j]);
 }

static __inline__ void
_27022015_upp_last_col(nmod_mat_t A,slong last_col,mp_limb_t** e)
 {
  mp_limb_t m=e[last_col][last_col]-1;
  slong j;
  for(j=last_col;j--;)
   e[j][last_col] &= m;
 }

static __inline__ void
_27022015_upp_easy(nmod_mat_t A,mp_limb_t mask)
// todo: replace division with shift
 {
  slong last_i=A->r-1,i,j;
  mp_limb_t** e=A->rows;
  for(i=1;i<last_i;i++)
   {
    mp_limb_t* sP=e[i]+i;
    mp_limb_t s=*sP;
    slong v_len=A->r-i;
    if(1==s)
     for(j=i;j--;)
      {
       mp_limb_t* tP=e[j]+i;
       mp_limb_t t_ori=*tP;
       if(t_ori)
        {
         t_ori &= mask;
         *tP = t_ori;
         if(t_ori)
          tmod_vec_muladd( tP, sP, v_len, -t_ori );
        }
      }
    else
     for(j=i;j--;)
      {
       mp_limb_t* tP=e[j]+i;
       mp_limb_t t_ori=*tP;
       if(t_ori >= s)
        {
         t_ori &= mask;
         *tP = t_ori;
         if(t_ori>=s)
          {
           t_ori /= s;
           tmod_vec_muladd( tP, sP, v_len, -t_ori );
          }
        }
      }
   }
  _27022015_upp_last_col(A,last_i,e);
 }

static __inline__ void
_27022015_upp_hard(nmod_mat_t A)
// todo: replace division with shift
 {
  slong last_i=A->r-1,i,j;
  mp_limb_t** e=A->rows;
  for(i=1;i<last_i;i++)
   {
    mp_limb_t* sP=e[i]+i;
    mp_limb_t s=*sP;
    slong v_len=A->r-i;
    if(1==s)
     for(j=i;j--;)
      {
       mp_limb_t* tP=e[j]+i;
       mp_limb_t t_ori=*tP;
       if(t_ori)
        tmod_vec_muladd( tP, sP, v_len, -t_ori );
      }
    else
     for(j=i;j--;)
      {
       mp_limb_t* tP=e[j]+i;
       mp_limb_t t_ori=*tP;
       if(t_ori >= s)
        {
         t_ori /= s;
         tmod_vec_muladd( tP, sP, v_len, -t_ori );
        }
      }
   }
  _27022015_upp_last_col(A,last_i,e);
 }

static __inline__ void
_27022015_upp(fmpz_t d,nmod_mat_t A)
 {
  HNF_mod_t_det_no_zeros(d,A);
  if( fmpz_cmp_ui( d, (UWORD(1)<<(FLINT_BITS-1))-1 ) <= 0 )
   {
    mp_limb_t mask=fmpz_to_t(d)-1;
    if(mask)
     _27022015_upp_easy(A,mask);
   }
  else
   _27022015_upp_hard(A);
 }

int
nmod_mat_HNF_mod_t(nmod_mat_t A,fmpz_t d)
/*
Put A into row-echelon form modulo T. Under-diagonal part of result undefined.

If one or more diagonal entries of result equal T, or diagonal is all ones, 
 then only diagonal elements are defined, and 1 returned
 
Else diagonal and above of result are HNF modulo T, and 0 is returned
 
d on exit equals product of A diagonal entries 
*/
 {
  #if LOUD_ROW_OP
   gmp_printf("before nmod_mat_HNF_mod_t():\n");
   nmod_mat_print_pretty(A);
   gmp_printf("taken modulo 64:\n");
   _27022015_print_diag_and_above(A,-1,64);
  #endif
  slong i;
  int zero_found=0;
  mp_limb_t* scrtch=(mp_limb_t*)malloc( A->r * sizeof(mp_limb_t) );
  for(i=0;i<A->r;i++)
   {
    zero_found |= _27022015_diag_and_down(A,i,scrtch);
    #if LOUD_ROW_OP
     gmp_printf("After cleanin column %d\n",i);
     _27022015_print_diag_and_above(A,i,64);
    #endif
   }
  if(zero_found)
   {
    fmpz_set_ui(d,1);
    fmpz_mul_2exp(d, d, FLINT_BITS);
   }
  else
   _27022015_upp(d,A);
  free(scrtch);
  #if LOUD_ROW_OP
   flint_printf("nmod_mat_HNF_mod_t() product:"); fmpz_print(d);
   flint_printf("\nnmod_mat_HNF_mod_t() return code: %d\n",zero_found);
   flint_printf("after nmod_mat_HNF_mod_t():\n");
   _27022015_print_diag_and_above(A,A->r,0);
  #endif
  return zero_found;
 }

#undef EMPTY_CACHE_CELL
#undef SENIOR_BIT_MASK
#undef HEAVY_CACHE_SIZE

#undef ASSERT
