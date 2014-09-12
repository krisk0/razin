// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#undef NDEBUG
#include <assert.h>
#include <flint/flint.h>

void
fmpz_mat_det_hermitian_decomposition(fmpz_t r, const fmpz_mat_t a)
 {
  if(a->r < 5)
   fmpz_mat_det_cofactor(r, a);
  else
   {
    fmpz_mat_t d;
    fmpz_t w; fmpz_init(w);
    int free_d=fmpz_mat_hermitian_decomposition_2(d,w,a);
    if( d->r )
     {
      // successfully decomposed, det a = (det d) * w
      fmpz_mat_det_odd(r,d);
      if(free_d)
       {
        fmpz_mat_clear(d);
        fmpz_mul(r,r,w);
       }
     }
    else
     {
      // singular or evil matrice a, w=2**126
      mpz_t m; mpz_init(m);
      __mpz_struct* z=_fmpz_promote_val(w);
      fmpz_mat_det_suspected_zero(m,a,z);
      _fmpz_demote_val(w);
      fmpz_set_mpz(r,m);
      mpz_clear(m);
     }
    fmpz_clear(w);
   }
 }

#undef NDEBUG
