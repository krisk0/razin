// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "../ulong_extras/profile_.h"
#include <flint/flint.h>

int fmpz_mat_hermitian_decomposition_2_64(fmpz_mat_t b,fmpz_t r, 
  const fmpz_mat_t m);

void
fmpz_mat_det_hermitian_decomposition(fmpz_t r, const fmpz_mat_t a)
 {
  if(a->r < 5)
   fmpz_mat_det_cofactor(r, a);
  else
   {
    fmpz_mat_t d;
    fmpz_t w; fmpz_init(w);
    #if SHOW_TIME
     flint_printf("\n\n");
    #endif
    MARK_TIME(t_deco);
    int free_d=fmpz_mat_hermitian_decomposition_2_64(d,w,a);
    DUMP_TIME("_hermitian_decomposition_2()",t_deco);
    if( d->r )
     {
      // successfully decomposed, det a = (det d) * w
      MARK_TIME(t0);
      fmpz_mat_det_odd(r,d);
      if(free_d)
       {
        fmpz_mat_clear(d);
        fmpz_mul(r,r,w);
       }
      DUMP_TIME("_det_odd()",t0);
     }
    else
     {
      // singular or evil matrice a, w is a degree of 2
      MARK_TIME(t1);
      mpz_t m; mpz_init(m);
      __mpz_struct* z=_fmpz_promote_val(w);
      fmpz_mat_det_suspected_zero(m,a,z);
      _fmpz_demote_val(w);
      fmpz_set_mpz(r,m);
      mpz_clear(m);
      DUMP_TIME("_det_suspected_zero()",t1);
     }
    fmpz_clear(w);
   }
 }

#undef NDEBUG
