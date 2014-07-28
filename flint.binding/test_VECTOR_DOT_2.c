#include <assert.h>
#include <flint/flint.h>
#include "C/ulong_extras/longlong_.h"
#include "C/ulong_extras/ulong_extras_.h"
#include "C/ulong_extras/mulmod_preinv_4arg.c"
#include "C/nmod_mat/nmod_mat_.h"

#define VECTOR_DOT_2_slow( r, w,x,y,z, M )   \
 {                                            \
  VECTOR_DOT_HEAD(w,x);                        \
  VECTOR_DOT_BODY(y,z);                           \
  VECTOR_DOT_TAIL(r, M.n,M.ninv,M.norm );            \
 }

flint_rand_t st;
p_k_pk_t P;
nmod_t M;
mp_limb_t A[4];

mp_limb_t R()
 {
  return n_randlimb(st) % M.n;
 }

#define TEST_slave(r0,r1, w,x,y,z, M)   \
    VECTOR_DOT_2_slow(r0, w,x,y,z, M ); \
    VECTOR_DOT_2(r1, w,x,y,z, M );      \
    if(r0 != r1)      \
     {                \
      flint_printf("n=%wu\n",M.n);      \
      flint_printf("w=%wu x=%wu y=%wu z=%wu\n",w,x,y,z);  \
      flint_printf("good=%wu baad=%wu\n",r0,r1);          \
      assert(0);      \
     }

#define TEST(r0,r1, w,x,y,z, M)         \
 A[0]=w; A[1]=x;       \
 A[2]=y; A[3]=z;        \
 TEST_slave(r0,r1, A[0],A[1],A[2],A[3], M) 
 
#define h0 UWORD(2626333620865894092)
#define h1 UWORD(9888751871984654438)
#define h2 UWORD(10486523025935645524)
#define h3 UWORD(588054462113098396)

int main()
 {
  P.p=3;
  slong i;
  mp_limb_t r0,r1;
  
  flint_randinit(st);
  init__p_k_pk__and__nmod( &P, &M );
  
  TEST(r0,r1, 0,0,0,0, M);
  TEST(r0,r1, 0,0,0,1, M);
  TEST(r0,r1, 0,0,1,1, M);
  TEST(r0,r1, 0,1,1,1, M);
  TEST(r0,r1, 1,1,1,1, M);
  TEST(r0,r1, 10,11,12,13, M);
  TEST(r0,r1, h0,1,1,1, M);

  TEST(r0,r1, h0, 0, 0, 0, M)
  TEST(r0,r1, h0, 0, 0,h3, M)
  TEST(r0,r1, h0, 0,h2,h3, M)
  TEST(r0,r1, h0,h1,h2,h3, M)
  
  TEST(r0,r1, h0,1+h0,1,1, M);
  TEST(r0,r1, h0,1+h0,2+h0,3+h0, M);
  
  for(i=100;i--;)
   {
    A[0]=R();  A[1]=R();
    A[2]=R();  A[3]=R();
    TEST_slave( r0,r1, A[0],A[1],A[2],A[3], M )
   }
  flint_printf("Test passed\n");
 }
