#include <flint/nmod_mat.h>

void 
nmod_mat_init_square_2arg(nmod_mat_t mat, slong dim)
// initialize square matrice with random values, do not touch mod structure
 {
  slong i;
  mp_limb_t* e = mat->entries = flint_malloc(dim * dim * sizeof(mp_limb_t));
  mat->rows = flint_malloc(dim * sizeof(mp_limb_t *));
  for (i = 0; i < dim; i++, e += dim)
    mat->rows[i] = e;
  mat->r = dim;
  mat->c = dim;
 }
