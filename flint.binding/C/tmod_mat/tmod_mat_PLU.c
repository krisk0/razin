// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

// Subroutines in this file are machine-dependent and known to work on amd64
//  with GMP

static __inline__ mp_limb_t
tmod_symm_abs( mp_limb_t y )
// return y or -y, whichever is smaller
 {
  if( y & 0x8000000000000000 )
   return -y;
  return y;
 }

static __inline__ mp_limb_t 
tmod_mat_PLU_find_nonzero(tmod_mat_t S, long* P, long col)
/*
put minimal in abs value element into S[col,col] by swapping rows
on success return inverse of S[col,col]
on failure return 0
*/
 {
  long m=S->r, min_i=col;
  mp_limb_t** a=S->rows;
  while(min_i<m)
   {
    if( a[min_i][col] & 1 )
     break;
    min_i ++;
   }
  if(min_i == m)
   return 0;
  long i,min_s=tmod_symm_abs( a[min_i][col] );
  mp_limb_t c;
  for(i=min_i+1;i<m;i++)
   {
    c=a[i][col];
    if( c & 1 )
     {
      c = tmod_symm_abs(c);
      if(c < min_s)
       {
        min_s=c; min_i=i;
       }
     }
   }
  if(min_i != col)
   {
    MP_PTR_SWAP( a[min_i], a[col] );
    i=P[min_i]; P[min_i]=P[col]; P[col]=i;
   }
  if(min_s == a[col][col])
   return t_invmod(min_s);
  return -t_invmod(min_s);
 }

static __inline__ void 
tmod_vec_scalar_addmul( mp_limb_t* tgt, mp_limb_t* sou, long L, 
 mp_limb_t Q )
 {
  while(L--)
   *tgt++ += (*sou++) * Q;
 }

// TODO: use block recursion for big m
long 
tmod_mat_PLU_mod_machine_word(long* PR, tmod_mat_t S)
/*
S: matrice with m rows and m-1 columns over residue ring modulo 2**64

attempt to find PLU factorisation of S, such that

P*original S = L * U

lower line of permutation P stored into PR, followed by array R such that
R[i] * diagonal(U)[i] = 1 modulo 2**64

return 1 on success, 0 on failure

S modified. On success S contains LU in compressed form, just like A after FLINT
 nmod_mat_lu(P,A,rank_check)

This subroutine is machine-dependent. Known to work on amd64, don't know what
 will happen on other arch
*/
 {
  long m=S->r, n=S->c, i,row,row_plus_1,length;
  for(i=m;i--;)
   PR[i]=i;                // initialize permutation P
  mp_limb_t** a=S->rows;
  mp_limb_t* alpha,* betta;
  mp_limb_t d,S_i_row,e;
  for(row=0;row<n;row++)
   {
    d=tmod_mat_PLU_find_nonzero( S, PR, row );
    if( !d )
     return 0;    // failed to find invertible element in this column, goodbye
    PR[m+row]=(long)d; // not good if mp_limb_t is larger than long
    row_plus_1 = row + 1;
    length = n - row_plus_1;
    alpha = a[row] + row_plus_1;
    for(i=row_plus_1;i<m;i++)
     {
      betta = a[i] + row_plus_1;
      S_i_row=betta[-1];
      if( S_i_row )
       {
        e = S_i_row * d;
        if(length != 0)
         tmod_vec_scalar_addmul( betta, alpha, length, -e );
        betta[-1] = e;
       }
     }
   }
  return 1;
 }
