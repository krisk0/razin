// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

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
  // transpose
  slong i,j,k,m=LU->r;
  slong cc=m-1;
  mp_limb_t* on_tgt,* on_sou;
  for(i=cc;i--;)
   {
    on_tgt = R->entries + i;   // will form column no. i
    on_sou = LU->rows[i+1];   // upper-left entry of L is at 1st row of LU
    // transfer i+1 pieces of L
    on_tgt[0] = on_sou[0];
    for(j=i;j--;)
     {
      on_tgt += cc;
      on_sou += 1;
      on_tgt[0] = on_sou[0];
     }
    on_sou = LU->rows[i]+i;
    on_tgt += cc;
    k = m-i-2;
    // transfer k+1 elements of W to on_tgt, on_tgt+cc, ...
    on_tgt[0] = on_sou[0];
    for(j=k;j--;)
     {
      on_tgt += cc;
      on_sou += 1;
      on_tgt[0] = on_sou[0];
     }
   }
 }
