#include <flint/nmod_mat.h>

void 
nmod_mat_init_3arg(nmod_mat_t mat, slong r, slong c)
// initialize matrice with random values, do not touch mod structure
 {
  slong i;
  mp_limb_t* e = mat->entries = flint_malloc(r * c * sizeof(mp_limb_t));
  mat->rows = flint_malloc(r * sizeof(mp_limb_t *));
  for (i = 0; i < r; i++, e += c)
    mat->rows[i] = e;
  mat->r = r;
  mat->c = c;
 }
