// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <flint/flint.h>
#include "../ulong_extras/ulong_extras_.h"
#include "../ulong_extras/longlong_.h"
#include "nmod_mat_.h"

#define TINY_DOT 0

/*
clang miscompiles this subroutine: 4x4x100 multiplication takes 
1.64 times longer under clang
*/
void
nmod_mat_mul_pk_classical(nmod_mat_t R,nmod_mat_t A,nmod_mat_t B)
 {
  slong i_max=R->r,j_max=R->c,k_max=B->r,i,j,k;
  mp_limb_t** const Rrows=R->rows;
  mp_limb_t** const Arows=A->rows;
  mp_limb_t** const Brows=B->rows;
  mp_limb_t* rho;
  mp_limb_t* alpha;
  mp_limb_t* betta;
  const mp_limb_t n=R->mod.n;
  const mp_limb_t ninv=R->mod.ninv;
  #if SPEEDUP_NMOD_RED3 && !TINY_DOT
   const mp_limb_t two_128_mod_n=R->mod.norm;
  #endif
  for(i=0;i<i_max;i++)
   {
    rho=Rrows[i];
    alpha=Arows[i];
    for(j=0;j<j_max;j++)
     {
      #if TINY_DOT
       VECTOR_DOT_HEAD_tiny(alpha[0],Brows[0][j], n,ninv);
       for(k=1;k<k_max;k++)
        VECTOR_DOT_BODY_tiny(alpha[k],Brows[k][j], n,ninv);
       VECTOR_DOT_TAIL_tiny(rho[j]);
      #else
       VECTOR_DOT_HEAD_greedy( alpha[0], Brows[0][j] );
       for(k=1;k<k_max;k++)
        VECTOR_DOT_BODY_greedy( alpha[k], Brows[k][j] );
       #if SPEEDUP_NMOD_RED3
        VECTOR_DOT_TAIL( rho[j], n,ninv,two_128_mod_n );
       #else
                                                // for 4x4x100 multiplication,
        VECTOR_DOT_TAIL_3arg( rho[j], n,ninv ); //  17% slower
       #endif
      #endif
     }
   }
 }

/*
asm code for this subroutine (made with gcc-4.7.3 -O3 -march=native) contains 
 line that look like garbage to me: xchg   %ax,%ax
0000000000400b20 <nmod_mat_mul_pk_classical>:
  400b20:	41 57                	push   %r15
  400b22:	41 56                	push   %r14
  400b24:	41 55                	push   %r13
  400b26:	41 54                	push   %r12
  400b28:	55                   	push   %rbp
  400b29:	53                   	push   %rbx
  400b2a:	48 8b 47 08          	mov    0x8(%rdi),%rax
  400b2e:	48 c7 44 24 d8 00 00 	movq   $0x0,-0x28(%rsp)
  400b35:	00 00 
  400b37:	48 8b 76 18          	mov    0x18(%rsi),%rsi
  400b3b:	48 8b 6f 20          	mov    0x20(%rdi),%rbp
  400b3f:	4c 8b 7f 28          	mov    0x28(%rdi),%r15
  400b43:	48 89 44 24 e0       	mov    %rax,-0x20(%rsp)
  400b48:	48 8b 47 10          	mov    0x10(%rdi),%rax
  400b4c:	4c 8b 5a 08          	mov    0x8(%rdx),%r11
  400b50:	48 89 74 24 f8       	mov    %rsi,-0x8(%rsp)
  400b55:	48 8b 5a 18          	mov    0x18(%rdx),%rbx
  400b59:	48 89 44 24 e8       	mov    %rax,-0x18(%rsp)
  400b5e:	48 8b 47 18          	mov    0x18(%rdi),%rax
  400b62:	4c 8b 74 24 e8       	mov    -0x18(%rsp),%r14
  400b67:	48 8b 7f 30          	mov    0x30(%rdi),%rdi
  400b6b:	48 89 44 24 f0       	mov    %rax,-0x10(%rsp)
  400b70:	49 c1 e6 03          	shl    $0x3,%r14
  400b74:	48 83 7c 24 e0 00    	cmpq   $0x0,-0x20(%rsp)
  400b7a:	48 89 7c 24 d0       	mov    %rdi,-0x30(%rsp)
  400b7f:	0f 8e fd 00 00 00    	jle    400c82 <nmod_mat_mul_pk_classical+0x162>
  400b85:	48 8b 54 24 d8       	mov    -0x28(%rsp),%rdx
  400b8a:	48 8b 44 24 f0       	mov    -0x10(%rsp),%rax
  400b8f:	48 83 7c 24 e8 00    	cmpq   $0x0,-0x18(%rsp)
  400b95:	4c 8b 2c d0          	mov    (%rax,%rdx,8),%r13
  400b99:	48 8b 44 24 f8       	mov    -0x8(%rsp),%rax
  400b9e:	4c 8b 14 d0          	mov    (%rax,%rdx,8),%r10
  400ba2:	0f 8e c4 00 00 00    	jle    400c6c <nmod_mat_mul_pk_classical+0x14c>
  400ba8:	4c 8b 23             	mov    (%rbx),%r12
  400bab:	45 31 c9             	xor    %r9d,%r9d
  400bae:	66 90                	xchg   %ax,%ax
  400bb0:	49 8b 02             	mov    (%r10),%rax
  400bb3:	4b f7 24 0c          	mulq   (%r12,%r9,1)
  400bb7:	49 83 fb 01          	cmp    $0x1,%r11
  400bbb:	49 89 c0             	mov    %rax,%r8
  400bbe:	48 89 d7             	mov    %rdx,%rdi
  400bc1:	0f 8e c6 00 00 00    	jle    400c8d <nmod_mat_mul_pk_classical+0x16d>
  400bc7:	31 f6                	xor    %esi,%esi
  400bc9:	b9 01 00 00 00       	mov    $0x1,%ecx
  400bce:	66 90                	xchg   %ax,%ax
  400bd0:	48 8b 14 cb          	mov    (%rbx,%rcx,8),%rdx
  400bd4:	49 8b 04 ca          	mov    (%r10,%rcx,8),%rax
  400bd8:	48 83 c1 01          	add    $0x1,%rcx
  400bdc:	4a f7 24 0a          	mulq   (%rdx,%r9,1)
  400be0:	49 01 c0             	add    %rax,%r8
  400be3:	48 11 d7             	adc    %rdx,%rdi
  400be6:	48 83 d6 00          	adc    $0x0,%rsi
  400bea:	4c 39 d9             	cmp    %r11,%rcx
  400bed:	75 e1                	jne    400bd0 <nmod_mat_mul_pk_classical+0xb0>
  400bef:	48 83 fe 01          	cmp    $0x1,%rsi
  400bf3:	76 14                	jbe    400c09 <nmod_mat_mul_pk_classical+0xe9>
  400bf5:	48 8b 44 24 d0       	mov    -0x30(%rsp),%rax
  400bfa:	48 f7 e6             	mul    %rsi
  400bfd:	31 f6                	xor    %esi,%esi
  400bff:	49 01 c0             	add    %rax,%r8
  400c02:	48 11 d7             	adc    %rdx,%rdi
  400c05:	48 83 d6 00          	adc    $0x0,%rsi
  400c09:	48 85 f6             	test   %rsi,%rsi
  400c0c:	74 10                	je     400c1e <nmod_mat_mul_pk_classical+0xfe>
  400c0e:	48 89 f8             	mov    %rdi,%rax
  400c11:	48 29 e8             	sub    %rbp,%rax
  400c14:	48 39 ef             	cmp    %rbp,%rdi
  400c17:	48 0f 47 f8          	cmova  %rax,%rdi
  400c1b:	48 29 ef             	sub    %rbp,%rdi
  400c1e:	48 89 f8             	mov    %rdi,%rax
  400c21:	48 29 e8             	sub    %rbp,%rax
  400c24:	48 39 ef             	cmp    %rbp,%rdi
  400c27:	48 0f 47 f8          	cmova  %rax,%rdi
  400c2b:	4c 89 f8             	mov    %r15,%rax
  400c2e:	48 f7 e7             	mul    %rdi
  400c31:	4c 01 c0             	add    %r8,%rax
  400c34:	48 11 fa             	adc    %rdi,%rdx
  400c37:	48 83 c2 01          	add    $0x1,%rdx
  400c3b:	48 0f af d5          	imul   %rbp,%rdx
  400c3f:	49 29 d0             	sub    %rdx,%r8
  400c42:	49 8d 14 28          	lea    (%r8,%rbp,1),%rdx
  400c46:	4c 39 c0             	cmp    %r8,%rax
  400c49:	4c 0f 46 c2          	cmovbe %rdx,%r8
  400c4d:	4c 89 c0             	mov    %r8,%rax
  400c50:	48 29 e8             	sub    %rbp,%rax
  400c53:	49 39 e8             	cmp    %rbp,%r8
  400c56:	4c 0f 43 c0          	cmovae %rax,%r8
  400c5a:	4f 89 44 0d 00       	mov    %r8,0x0(%r13,%r9,1)
  400c5f:	49 83 c1 08          	add    $0x8,%r9
  400c63:	4d 39 f1             	cmp    %r14,%r9
  400c66:	0f 85 44 ff ff ff    	jne    400bb0 <nmod_mat_mul_pk_classical+0x90>
  400c6c:	48 83 44 24 d8 01    	addq   $0x1,-0x28(%rsp)
  400c72:	48 8b 54 24 e0       	mov    -0x20(%rsp),%rdx
  400c77:	48 39 54 24 d8       	cmp    %rdx,-0x28(%rsp)
  400c7c:	0f 85 03 ff ff ff    	jne    400b85 <nmod_mat_mul_pk_classical+0x65>
  400c82:	5b                   	pop    %rbx
  400c83:	5d                   	pop    %rbp
  400c84:	41 5c                	pop    %r12
  400c86:	41 5d                	pop    %r13
  400c88:	41 5e                	pop    %r14
  400c8a:	41 5f                	pop    %r15
  400c8c:	c3                   	retq   
  400c8d:	49 89 c0             	mov    %rax,%r8
  400c90:	48 89 d7             	mov    %rdx,%rdi
  400c93:	eb 89                	jmp    400c1e <nmod_mat_mul_pk_classical+0xfe>
  400c95:	66 66 2e 0f 1f 84 00 	data32 nopw %cs:0x0(%rax,%rax,1)
  400c9c:	00 00 00 00 
*/
