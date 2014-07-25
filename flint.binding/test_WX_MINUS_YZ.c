#include <assert.h>
#include <flint/flint.h>
#include "C/ulong_extras/longlong_.h"
#include "C/ulong_extras/ulong_extras_.h"
#include "C/ulong_extras/mulmod_preinv_4arg.c"

#define WX_MINUS_YZ_slow( r, w,x,y,z, n,i ) \
 r=n_submod( n_mulmod_preinv_4arg(w,x,n,i), n_mulmod_preinv_4arg(y,z,n,i), n );

flint_rand_t st;
p_k_pk_t P;
nmod_t M;

mp_limb_t R()
 {
  return n_randlimb(st) % M.n;
 }

int main()
 {
  P.p=3;
  slong i;
  mp_limb_t w,x,y,z,r0,r1;
  
  flint_randinit(st);
  init__p_k_pk__and__nmod( &P, &M );
  
  for(i=100;i--;)
   {
    w=R();  x=R();
    y=R();  z=R();
    WX_MINUS_YZ_slow(r0, w,x,y,z, M.n,M.ninv );
    WX_MINUS_YZ(r1, w,x,y,z, M.n,M.ninv,M.norm );
    if(r0 != r1)
     {
      flint_printf("n=%wu\n",M.n);
      flint_printf("w=%wu x=%wu y=%wu z=%wu\n",w,x,y,z);
      flint_printf("good=%wu baad=%wu\n",r0,r1);
      assert(0);
     }
   }
  flint_printf("Test passed\n");
 }
