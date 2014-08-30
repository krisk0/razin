// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>

void tmod_mat_mul_diag_tril( tmod_mat_t A, const mp_limb_t* D )
// A := D*A, A square lower-triangular; D diagonal, represented by array
 {
  slong m=A->r,i,j;
  mp_limb_t d;
  mp_limb_t* e;
  mp_limb_t** a=A->rows;
  for(i=0;i<m;i++) // access memory in ascending order if rows are unshifted
   {
    d=D[i];
    if(d!=1)
     {
      e=a[i];
      for(j=i+1;j--;e++)
       *e *= d;    // access memory in ascending order
     }
   }
 }

void tmod_mat_mul_trilV_diag( tmod_mat_t A, const mp_limb_t* D )
// A := A*D, A square lower-triangular virgin; D diagonal, represented by array
 {
  slong m=A->r, i,j,k;
  slong m_plus=m+1;
  mp_limb_t d;
  mp_limb_t* e;
  mp_limb_t* a=A->entries;
  for(i=0;i<m;i++)
   {
    d=D[i];
    if(d!=1)
     {
      e=a + i*m_plus;    // point to diagonal
      k=m-i;            // column (generally) contains this many non-zeros
      for(j=k;j--;e += m)
       *e *= d;
     }
   }
 }

void 
tmod_mat_invert_W_part( tmod_mat_t R, const mp_limb_t* D )
/*
R on entry: m rows, m-1 cols, lower part W below diagonal
 
R virgin
 
R on exit: W part inverted, L part untouched, rows unshifted
 
          &  @                                        &   @
input R = 1  ~  W = 1    D = (1, 1/3)     output R =  1   ~
         -6  3     -6 3                               2 1/3
*/
 {
  tmod_mat_window_unsh_t C;
  slong m=R->r, i,j,k;
  slong cc=m-1;
  tmod_mat_t E;
  tmod_mat_init( E, cc, cc );
  mp_limb_t* Erow=E->entries;
  mp_limb_t* Wrow=R->rows[1];  // top-left entry of W is at row 1
  // un-compress W to E
  for(i=0;i<cc;i++,Erow += cc, Wrow += cc)
   {
    k=i+1;    // for i=cc-1 transfer all cc entries
    for(j=0;j<k;j++)
     Erow[j] = Wrow[j];
    //assert(Erow[i] * D[i] == 1 );
   }
  tmod_mat_mul_diag_tril( E, D );
  tmod_mat_to_window_unsh( C, E, 0, 0, cc, cc );
  tmod_mat_invert_tril( C );
  tmod_mat_window_unsh_clear( C );
  tmod_mat_mul_trilV_diag( E, D );
  // compress E to W
  Wrow=R->rows[1];
  Erow=E->entries;
  for(i=0;i<cc;i++, Erow += cc, Wrow += cc)
   {
    k=i+1;
    for(j=0;j<k;j++)
     Wrow[j] = Erow[j];
   }
  tmod_mat_clear( E );
 }

void
tmod_mat_invert_L_part( tmod_mat_t R, const tmod_mat_t S)
/*
S: m rows, m-1 cols, lower part L below diagonal
 
R: target, of same dimension, virgin
 
R on exit: L part of S inverted and placed into upper part; still virgin
 
          &  @      1              1                        -4 15
input S = 4  ~  L = 4  1     L.I= -4  1          output R =  ?  2
          7 -2      7 -2  1       15  2  1                   ?  ?
*/
 {
  tmod_mat_window_unsh_t C;
  slong m=R->r, i,j;
  slong cc=m-1;
  tmod_mat_t E;
  tmod_mat_init( E, m, m );
  mp_limb_t* Erow=E->entries;
  mp_limb_t* Lrow; 
  // un-compress L to E
  Erow[0] = 1;
  Erow += m;
  for (i=1; i<m; i++, Erow += m)
   {
    Lrow=S->rows[i];       // top-left entry of L is at row 1
    for (j=0;j<i;j++)     // row no. i wants i general elements, 
     Erow[j] = Lrow[j];
    Erow[i] = 1;         //                                     followed by 1
   }
  tmod_mat_to_window_unsh( C, E, 0, 0, m, m );
  tmod_mat_invert_tril( C );
  tmod_mat_window_unsh_clear( C );
  // compress-transpose E
  for (i=0;i<cc;i++)
   {
    Lrow=R->entries+i;
    j=i+1;
    Erow=E->rows[j];
    // put j entries starting at Erow into Lrow, Lrow+cc, Lrow+2*cc
    for ( ;j--; Lrow += cc)
     *Lrow = *Erow++;
   }
  tmod_mat_clear( E );
 }

void
tmod_mat_solver_3arg( tmod_mat_t R, long* PD, const tmod_mat_t LU )
/*
 PD, LU: jammed PR, jammed LU as returned by tmod_mat_PLU(), with m rows
  and m-1 columns
 
 R: matrice with un-shifted rows, same dimensions as LU
 
 invert transposed L and transposed Uu where Uu = U without zero row
 
 put result into R, L'=L.T.I jammed into upper part,
  U'=Uu.T.I jammed into lower part, as in example below
  
 PD: array ( ?, ?, ?, 1, 0xAAAAAAAAAAAAAAAB )
  
      1 -6       1             1 -6
 LU = 4  3   L = 4  1      U =    3
      7 -2       7 -2  1
     
 Denote with 1/3 number 3**(-1) modulo 2**64=0xAAAAAAAAAAAAAAAB

     -4 -15
 R =  1   2     - result for such input
      2 1/3
*/
 {
  // invert L while it is in place
  tmod_mat_invert_L_part( R, LU );
  // put upper part of LU into lower part of R
  slong i,j,k,m=LU->r;
  slong cc=m-1;
  mp_limb_t* on_tgt,* on_sou;
  for ( i=0; i<cc; i++ )
   {
    k = m-i-2;
    on_tgt = R->rows[1+i]+i;
    on_sou = LU->rows[i]+i;
    // transfer k+1 elements of W to on_tgt, on_tgt+cc, ...
    on_tgt[0] = on_sou[0];
    for ( ;k--; )
     {
      on_tgt += cc;
      on_sou += 1;
      on_tgt[0] = on_sou[0];
     }
   }
  // invert lower W part
  tmod_mat_invert_W_part( R, (mp_limb_t*) (PD+m) );
 }
