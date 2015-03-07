#if WANT_ASSERT_IN_HNF_NONSQ
 #include <assert.h>
#else
 #define assert(x) 
#endif

static __inline__ slong
DKryskov_fix_diagonal_tail(nmod_mat_t A,slong i,mp_limb_t n)
//return 1 iff n is not a multiple of diagonal tail
 {
  mp_limb_t* eP=&nmod_mat_entry(A,i,i);
  mp_limb_t e=*eP % n;
  if(e)
   {
    *eP = n_gcd(n,e);
    return 1;
   }
  *eP = n;
  return 0;
 }

static __inline__ mp_limb_t
DKryskov_last_col_episode( mp_limb_t pivot, mp_limb_t* other )
 {
  mp_limb_t x = *other % pivot;
  if(x)
   {
    pivot = n_gcd( pivot, x );
    *other = 0;
   }
  return pivot;
 }

static __inline__ void
DKryskov_HNF_zap_below_last_col(nmod_mat_t A,slong col,slong hint,mp_limb_t n)
 {
  slong i,m;
  if(hint == -3)
   {
    for(i=A->r-1; i>col; i--)
     nmod_mat_entry(A,i,col) = 0;
   }
  if(hint < 0)
   return;
  mp_limb_t mu = nmod_mat_entry(A,col,col);
  mp_limb_t nu = n % mu;
  if(nu)
   nu = n_gcd( n, nu );
  else
   nu = mu;
  // hint = line to be zapped 1st 
  nu = DKryskov_last_col_episode( nu, &nmod_mat_entry(A,hint,col) );
  i=A->r-1;
  while( (nu > 1) && (i>col) )
   {
    nu = DKryskov_last_col_episode( nu, &nmod_mat_entry(A,i,col) );
    i--;
   }
  // nu=1 or i=col
  for( ; i>col; i--)
   nmod_mat_entry(A,i,col) = 0;
  if(mu != nu)
   nmod_mat_entry(A,col,col) = nu;
 }

static __inline__ void
DKryskov_nmod_1_lower_nonsq(nmod_mat_t A,slong col,slong j,mp_limb_t n)
/*
A[col,col] is known to be 1
zap column col starting from row j
operate modulo n
*/
 {
    #if BUG0_nmod_mat_HNF
     printf("DKryskov_nmod_1_lower() start col=%ld j=%ld n=%lu\n",col,j,n);
     nmod_mat_print(A);
    #endif
  slong m=A->r;
  mp_limb_t* sP=A->rows[col]+col;
  assert( 1 == *sP );
  //TODO: skip counting elements of column col, to save time
  slong v_len = A->c-col;
  while(j < m)
   {
    mp_limb_t* tP=A->rows[j]+col;
    mp_limb_t t_ori=*tP;
    if(t_ori)
     {
      mp_limb_t t_upd=t_ori % n;
      if(t_upd)
       {
        #if BUG0_nmod_mat_HNF
         printf("q=%lu old vector=",n-t_upd);
         vec_print( "", tP, v_len );
        #endif
        _nmod_vec_scalar_addmul_nmod( tP, sP, v_len, n-t_upd, A->mod );
        #if BUG0_nmod_mat_HNF
         vec_print( "new vector=", tP, v_len );
         printf("\n");
        #endif
       }
      else
       *tP = 0; // TODO: this line can be deleted?
     }
    ++j;
   }
    #if BUG0_nmod_mat_HNF
     printf("DKryskov_nmod_1_lower() end\n");
     nmod_mat_print(A);
    #endif
 }

static __inline__ void
DKryskov_HNF_zap_below(nmod_mat_t A,slong col,slong hint,mp_limb_t n,
  mp_limb_t* scrtch)
 {
  if(hint == -2)
   {
    nmod_mat_entry(A,col,col) = n;
    return;
   }
  if(hint == -3)
   DKryskov_nmod_1_lower_nonsq(A,col,col+1,n);
  if(hint < 0)
   return;
  slong m = A->r, j;
  if(DKryskov_nmod_zero_line(A,col,hint,n,scrtch))
   {
    j = col+1;
    if(j < m)
     DKryskov_nmod_1_lower_nonsq(A,col,j,n);
    return;
   }
  j = col;
  while( ++j < m )
   {
    if( nmod_mat_entry(A,j,col) % n )
     {
      if(DKryskov_nmod_zero_line(A,col,j,n,scrtch))
       {
        if( ++j < m )
         DKryskov_nmod_1_lower_nonsq(A,col,j,n);
        return;
       }
     }
   }
 }

#if !WANT_ASSERT_IN_HNF_NONSQ
 #undef assert
#endif
