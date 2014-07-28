// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <assert.h>
#include <flint/flint.h>
#include "C/ulong_extras/longlong_.h"
#include "C/ulong_extras/ulong_extras_.h"
#include "C/ulong_extras/mulmod_preinv_4arg.c"
#include "C/nmod_mat/nmod_mat_.h"

#if 0
 #define PRINTF flint_printf
#else
 #define PRINTF
#endif

#define VECTOR_DOT_2_slow( r, w,x,y,z, M )   \
 {                                            \
  PRINTF("%wu * %wu + %wu * %wu\n", w,x,y,z ); \
  VECTOR_DOT_HEAD(w,x);                        \
  PRINTF("w*x=%wu %wu\n", Vmi,Vlo ); \
  VECTOR_DOT_BODY(y,z);                           \
  PRINTF("w*x + y*z=%wu %wu %wu\n",Vhi,Vmi,Vlo ); \
  VECTOR_DOT_TAIL(r, M.n,M.ninv,M.norm );            \
 }

 #define minus1 UWORD(-1)
 
 #define VECTOR_DOT_2_diag(tgt, x0,y0, x1,y1, mod)   \
  {                                                    \
   register mp_limb_t V2_dx asm ("rdx");                    \
   register mp_limb_t V2_ax asm ("rax");                       \
   mp_limb_t V2_me,V2_lo; /* 1 register for n + theese 2 = 3 */  \
   register mp_limb_t n=mod.n;                                    \
   mp_limb_t d_Vhi=minus1,d_sub0=minus1,d_sub1=minus1;                   \
   asm                                                            \
    (                                                             \
        "movq               $0,%q4\n\t" \
        "mov  %q7,%q1\n\t"                                        \
        "mulq %q8\n\t"                                            \
        "movq  %q1,%q3\n\t"                                        \
        "mov  %q0,%q2\n\t"                                        \
        "mov  %q9,%q1\n\t"                                        \
        "mulq %q10\n\t"                                            \
        "add  %q1,%q3\n\t"                                        \
        "adc  %q0,%q2\n\t"                                        \
        "jnc  1f\n\t"                                             \
        /* maybe subtract n from V2_me */                         \
        "movq               $1,%q4\n\t" \
        "xor  %q1,%q1\n\t"                                        \
        "cmp  %q2,%q11\n\t"                                        \
        "cmovb %q11,%q1\n\t"   /* if V2_me >= n then ax=n */      \
        "movq              %q1,%q5\n\t" \
        "sub  %q1,%q2\n\t"                                        \
        /* always subtract n from V2_me */                        \
        "sub  %q11,%q2\n\t"                                        \
     "1: xor  %q0,%q0\n\t"                                        \
        "cmp  %q2,%q11\n\t"                                        \
        "cmovb %q11,%q0\n\t"   /* if V2_me >= n then dx=n */      \
        "sub %q0,%q2\n\t"                                         \
        "movq              %q2,%q6\n\t" \
     : "=&d" (V2_dx), "=&a" (V2_ax), "=&r" (V2_me), "=&r" (V2_lo), \
       "+m" (d_Vhi), "+m" (d_sub0), "+m" (d_sub1)   \
     : "m" (x0),"m" (y0), "m" (x1),"m" (y1), "r" (n)            \
     ); \
   flint_printf("Vhi=%wu sub0=%wu tmp1=%wu\n",d_Vhi,d_sub0,d_sub1); \
   flint_printf("2 limb~ %wu %wu\n",V2_me,V2_lo); \
   NMOD_RED2_pk_4arg(V2_me,V2_lo, n,mod.ninv);               \
   tgt=V2_lo;                                             \
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
      flint_printf("good=%wu baad=%wu\n n=%wu p=%wu k=%wu\n",r0,r1,M.n, \
       P.p,P.k); \
      assert(0);      \
     }

#define TEST(r0,r1, w,x,y,z, M)         \
 A[0]=w; A[1]=x;       \
 A[2]=y; A[3]=z;        \
 TEST_slave(r0,r1, A[0],A[1],A[2],A[3], M) 

#define TEST_careful(r0,r1, w,x,y,z, M)        \
 A[0]=UWORD(w) % M.n; A[1]=UWORD(x) % M.n;   \
 A[2]=UWORD(y) % M.n; A[3]=UWORD(z) % M.n;  \
 TEST_slave(r0,r1, A[0],A[1],A[2],A[3], M); 
 
#define h0 UWORD(2626333620865894092)
#define h1 UWORD(9888751871984654438)
#define h2 UWORD(10486523025935645524)
#define h3 UWORD(588054462113098396)

void 
do_test()
 {
  slong i;
  mp_limb_t r0,r1;
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
  
  for(i=1000;i--;)
   {
    A[0]=R();  A[1]=R();
    A[2]=R();  A[3]=R();
    TEST_slave( r0,r1, A[0],A[1],A[2],A[3], M )
   }

  TEST_careful( r0,r1,
   13835058055282163708,13835058055282163710,
   13835058055282163705,13835058055282163710,
   M );

  TEST_careful( r0,r1,
   13844340530321392734,13993208285623501322,
   10359436948649161513,13802576360186490278,
   M );
 }

int main()
 {
  mp_limb_t z;

  flint_randinit(st);

  P.p=3;
  init__p_k_pk__and__nmod( &P, &M );
  do_test();

  P.p=5;
  init__p_k_pk__and__nmod( &P, &M );
  do_test();

  P.k=1;
  P.p_deg_k=P.p;
  count_leading_zeros_opt( z, P.p_deg_k );
  M.n=P.p<<z;
  invert_limb(M.ninv,M.n);
  z = - M.n;
  M.norm = n_mulmod_preinv_4arg( z,z, M.n,M.ninv );
  
  do_test();
  
  flint_printf("Test passed\n");
 }
