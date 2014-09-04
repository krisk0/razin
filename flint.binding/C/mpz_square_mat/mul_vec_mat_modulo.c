// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "mpz_square_mat_.h"

void mpz_square_mat_mul_vec_mat_modulo(mpz_ptr t,
  mpz_ptr v,const mpz_square_mat_t A,mpz_t m)
// t := v*A mod m
 {
  const slong n=A->r;
  slong i,j;
  mpz_ptr* ro=A->rows;
  mpz_ptr vP,aP;
  mpz_ptr tP;
  for(i=0,tP=t;i<n;i++,tP++)
   {
    mpz_set_ui(tP,0);
    for(vP=v,j=0;j<n;j++,vP++)
     {
      aP = ro[i]+j;
      mpz_addmul(tP, vP, aP);
     }
    mpz_tdiv_r(tP, tP, m);
   }
 }
