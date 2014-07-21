#include <flint/flint.h>
#include <flint/nmod_mat.h>

void 
nmod_mat_randfill(nmod_mat_t m, flint_rand_t f)
// fill matrice with random values, uniform distribution not promised
 {
  slong i;
  mp_limb_t n=m->mod.n;
  mp_limb_t* e=m->entries;
  for(i=m->r*m->c; i--; )
   *e++ = n_randbits( f, 64 ) % n;
 }
