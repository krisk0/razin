static __inline__ void
DKryskov_fix_diagonal_tail(nmod_mat_t A,slong i,mp_limb_t n)
 {
  mp_limb_t* eP=&nmod_mat_entry(A,i,i);
  mp_limb_t e=*eP % n;
  if(e)
   *eP = n_gcd(n,e);
  else
   *eP = n;
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
DKryskov_HNF_zap_below(nmod_mat_t A,slong col,slong hint,mp_limb_t n,
  mp_limb_t* scrtch)
 {
  if(hint == -2)
   {
    nmod_mat_entry(A,col,col) = n;
    return;
   }
  if(hint == -3)
   DKryskov_nmod_1_lower(A,col,col+1,n);
  if(hint < 0)
   return;
  slong m = A->r, j;
  if(DKryskov_nmod_zero_line(A,col,hint,n,scrtch))
   {
   }
  while( ++j < m )
   {
    
   }
 }

static __inline__ void
DKryskov_Gauss_upper(nmod_mat_t A)
 {
  assert(0);
 }
