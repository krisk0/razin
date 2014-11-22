// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "../ulong_extras/ulong_extras_.h"

void
tmod_mat_print_hex(char* m,const tmod_mat_t S)
 {
  slong i,j=FLINT_BITS;
  for(i=S->r*S->r;i--;)
   {
    mp_limb_t c=S->entries[i];
    if(c)
     {
      mp_limb_t count;
      count_leading_zeros_opt(count, c);
      if(count<j)
       j=count;
     }
   }
  // FLINT_BITS-j = max bit-length of S entries
  j=(FLINT_BITS-j+15)>>2;
  j += (!j);
  char fmt[10];
  sprintf(fmt,"%c%dMX ",'%',j);
  if(m)
   gmp_printf("%s\n",m);
  for(i=0;i<S->r;i++)
   {
    mp_limb_t* p=S->rows[i];
    for(j=0;j<S->r;j++)
     gmp_printf(fmt,p[j]);
    gmp_printf("\n");
   }
 }

static __inline__ void
_20141102_Lo(tmod_mat_t R, const tmod_mat_t S, slong m)
/*
 invert/transpose lower part of square S, put result into upper part of R.
 Destroy diagonal and lower part of R

 R virgin
*/ 
 {
  tmod_mat_window_unsh_t C;
  slong i,j,k;
  mp_limb_t* Erow=R->entries;
  mp_limb_t* Lrow;
  // fill lower part and diagonal of R
  Erow[0] = 1;
  Erow += m;
  for(i=1; i<m; i++, Erow += m)
   {
    Lrow=S->rows[i];       // top-left entry of L is at row 1
    for (j=0;j<i;j++)     // row no. i wants i general elements, 
     Erow[j] = Lrow[j];
    Erow[i] = 1;         //                                     followed by 1
   }
  //tmod_mat_print_hex("after moving lower part:",R);
  tmod_mat_to_window_unsh( C, R, 0, 0, m, m );
  tmod_mat_invert_tril( C );
  tmod_mat_window_unsh_clear( C );
  //tmod_mat_print_hex("after inverting lower part:",R);
  // transfer result to upper part of R
  for(i=1;i<m;i++)
   {
    Erow=R->entries+i;
    Lrow=R->rows[i];
    for(j=i;j--;Erow += m,Lrow++) // transfer i entries from row i to col i
     Erow[0] = Lrow[0];
   }
 }

#if 0
static __inline__ void
tmod_mat_mul_diag(tmod_mat_t A, const mp_limb_t* D)
// A := D*A, A square virgin; D diagonal, represented by array
 {
  slong m=A->r,i,j;
  mp_limb_t d;
  mp_limb_t* e;
  for(i=0;i<m;i++)
   {
    d=D[i];
    e=A->rows[i];
    for(j=m;j--;e++)
     e[0] *= d;    // access memory in ascending order
   }
 }
#endif

static __inline__ void
_20141102_Up(tmod_mat_t R, const tmod_mat_t S, mp_limb_t* dI, slong n)
/*
 invert/transpose upper part of square S, put result into diagonal and 
  lower part of R

 R virgin  
*/
 {
  {
   // transpose, multiply by Diag(dI)
   slong i,k;
   mp_limb_t* on_tgt,* on_sou;
   mp_limb_t Di;
   for(i=0; i<n; i++)
    {
     k = n-i-1;
     Di=dI[i];
     on_tgt = R->rows[i]+i;
     on_sou = S->rows[i]+i;
     // transfer k+1 elements of W to on_tgt, on_tgt+cc, ...
     on_tgt[0] = 1; // Di*on_sou[0];
     for( ;k--; )
      {
       on_tgt += n;
       on_sou += 1;
       on_tgt[0] = Di*on_sou[0];
      }
    }
  }
  /*
    t
   U  = W * D
   
    t '
   U    = D' * W'
  */
  // invert W now stored in lower part of R
  {
   tmod_mat_window_unsh_t C;
   tmod_mat_to_window_unsh(C, R, 0,0, n,n);
   tmod_mat_invert_tril(C);
   tmod_mat_window_unsh_clear(C);
  }
  tmod_mat_mul_diag_tril(R, dI);
 }

static __inline__ void
_20141102_unLU(tmod_mat_t R, const tmod_mat_t S, slong n)
/*
get back from LU representation

   a   b
S=
   c   d


   1   b       a   0
R=         *  
   0   1       c   d
*/
 {
  mp_limb_t rez;
  mp_limb_t* rho=R->entries;
  mp_limb_t* cOl;
  mp_limb_t* rOw;
  slong i,j,k;
  for(i=0;i<n;i++)
   {
    /*
    j<=i
    ( 1 s[i,i+1] s[i,i+2] ... s[i,n-1] ) * ( s[i,j] s[i+1,j] ... s[n-1,j] )
    */
    for(j=0;j<=i;j++)
     {
      cOl=S->rows[i]+j;
      rez=cOl[0];
      rOw=S->rows[i]+i;
      for(k=i+1;k<n;k++)
       {
        cOl += n;
        rOw ++;
        rez += rOw[0]*cOl[0];
       }
      *rho++ = rez;
     }
    /*
    j>i
    ( s[i,j] s[i,j+1] ... s[i,n-1] ) * ( s[j,j] s[j+1,j] ... s[n-1,j] )
    */
    for(j=i+1;j<n;j++)
     {
      rOw=S->rows[i]+j;
      cOl=S->rows[j]+j;
      rez=0;
      for(k=j;k<n;k++,rOw++,cOl+=n)
       rez += rOw[0]*cOl[0];
      *rho++ = rez;
     }
   }
 }

static __inline__ void
_20141102_shift_rows(tmod_mat_t R, tmod_mat_t Q, mp_limb_t* P, slong n, 
  int P_parity)
// multiply R by permutation matrice on the left. Use Q as scratch
 {
  // if P_parity is negative then P is not identity
  if( (1==P_parity) && is_identity_permutation(P,n) )
   return;
  {
   mp_limb_t** C=Q->rows;
   mp_limb_t** B=R->rows;
   memcpy(C,B,n*sizeof(mp_limb_t**));
   slong i;
   for(i=0;i<n;i++)
    B[ P[i] ]=C[ i ];
   // C destroyed
   tmod_mat_virginize(Q);
   // C restored to virgin state
   slong size=n*sizeof(mp_limb_t);
   for(i=0;i<n;i++)
    memcpy(C[i],B[i],size);
   // matrix contents copied
   R->rows=C; Q->rows=B;
  }
  {
   MP_PTR_SWAP(R->entries,Q->entries);
  }
 }

mp_limb_t
tmod_mat_invert_transpose(tmod_mat_t R, const tmod_mat_t S)
/*
S,R square, virgin

If S in inverible, R:=S inverted transposed in virgin state, return 1

Else return 0
*/
 {
  slong n=S->r;
  slong row_size=n*sizeof(mp_limb_t);
  memcpy(R->entries,S->entries,n*row_size);
  mp_limb_t* PR=flint_malloc(2*row_size);
  mp_limb_t d=0;
  if(tmod_mat_PLU_mod_machine_word(PR,R))
   {
    int p_parity=count_permutation_parity(PR,n);
    d=tmod_mat_diag_product_ZZ_ui(R)*p_parity;
    tmod_mat_t K; tmod_mat_init_fast(K,n,n);
    _20141102_Lo(K,R,n);
    _20141102_Up(K,R,PR+n,n);
    tmod_mat_virginize(R);
    _20141102_unLU(R,K,n);
    _20141102_shift_rows(R,K,PR,n,p_parity);
    tmod_mat_clear(K);
   }
  flint_free(PR);
  return d;
 }
