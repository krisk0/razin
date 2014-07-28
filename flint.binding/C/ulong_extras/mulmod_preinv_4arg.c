// This program is part of RAZIN
/******************************************************************************

    Copyright (C) 2009, 2010, 2012 William Hart

******************************************************************************/
// Copyright Денис Крыськов 2014

// This file contains a piece of FLINT mulmod_preinv.c owned by W.Hart

// Licence: GNU General Public License (GPL)

#include <flint/flint.h>

mp_limb_t
n_mulmod_preinv_4arg(mp_limb_t a, mp_limb_t b, mp_limb_t n, mp_limb_t ninv)
//castrated FLINT subroutine n_mulmod_preinv(): norm parameter removed
{
 #if MP_LIMB_BITS == 64 && defined (__amd64__)
   register mp_limb_t tgt asm ("rax");
   register mp_limb_t scr; // asm ("r8")
   asm 
    (
     "mov %q0,%q4\n\t"      // %q4 := n
     "mov %q2,%q1\n\t"
     "mulq %q3\n\t"         // a,b no longer needed
     "mov %q0,%q3\n\t"
     "mov %q1,%q2\n\t"
     "mov %q0inv,%q1\n\t" 
     "mulq %q0\n\t"         // ninv no longer needed
     "add %q2,%q1\n\t"
     "adc %q3,%q0\n\t"
     "add $1,%q0\n\t"       // gcc insists on using add rather than inc
     "imul %q4,%q0\n\t"
     "sub %q0,%q2\n\t"
     "lea (%q2,%q4,1),%q0\n\t"
     "cmp %q2,%q1\n\t"
     "cmovbe %q0,%q2\n\t"
     "mov %q2,%q1\n\t"      // tgt = %q2
     "sub %q4,%q1\n\t"      // tgt -= n
     "cmp %q4,%q2\n\t"      // ? n > %q2
     "cmovb %q2,%q1\n\t"    // if %q2 < n then tgt = %q2
     : "+d" (n), "=&a" (tgt), "+r" (a), "+r" (b), "=r" (scr)
     : "r" (ninv)
    );
   return tgt;
 #else
     mp_limb_t q0, q1, r, p;
 
     /* multiply */
     umul_ppmm(p, r, a, b);
     
     /* reduce mod n */
     umul_ppmm(q1, q0, ninv, p);
     add_ssaaaa(q1, q0, q1, q0, p, r);
 
     r -= (q1 + 1) * n;
 
     if (r >= q0)
      r += n;
 
     return (r < n ? r : r - n);
 #endif
}


/* asm code made with gcc 4.7.3 -O3:
   1be20:	49 89 d0             	mov    %rdx,%r8
   1be23:	48 89 f8             	mov    %rdi,%rax
   1be26:	48 f7 e6             	mul    %rsi
   1be29:	48 89 c7             	mov    %rax,%rdi
   1be2c:	48 89 d6             	mov    %rdx,%rsi
   1be2f:	48 89 c8             	mov    %rcx,%rax
   1be32:	48 f7 e2             	mul    %rdx
   1be35:	48 89 c1             	mov    %rax,%rcx
   1be38:	48 89 f8             	mov    %rdi,%rax
   1be3b:	48 01 f9             	add    %rdi,%rcx
   1be3e:	48 11 f2             	adc    %rsi,%rdx
   1be41:	48 83 c2 01          	add    $0x1,%rdx
   1be45:	49 0f af d0          	imul   %r8,%rdx
   1be49:	48 29 d0             	sub    %rdx,%rax
   1be4c:	4a 8d 14 00          	lea    (%rax,%r8,1),%rdx
   1be50:	48 39 c1             	cmp    %rax,%rcx
   1be53:	48 0f 46 c2          	cmovbe %rdx,%rax
   1be57:	48 89 c2             	mov    %rax,%rdx
   1be5a:	4c 29 c2             	sub    %r8,%rdx
   1be5d:	4c 39 c0             	cmp    %r8,%rax
   1be60:	48 0f 43 c2          	cmovae %rdx,%rax
   1be64:	c3                   	retq   

gcc 4.8.3 2 lines shorter:
   1bbe0:	49 89 d0             	mov    %rdx,%r8
   1bbe3:	48 89 f8             	mov    %rdi,%rax
   1bbe6:	48 f7 e6             	mul    %rsi
   1bbe9:	48 89 d6             	mov    %rdx,%rsi
   1bbec:	48 89 c7             	mov    %rax,%rdi
   1bbef:	48 89 c8             	mov    %rcx,%rax
   1bbf2:	48 f7 e2             	mul    %rdx
   1bbf5:	48 01 f8             	add    %rdi,%rax
   1bbf8:	48 11 f2             	adc    %rsi,%rdx
   1bbfb:	48 83 c2 01          	add    $0x1,%rdx
   1bbff:	49 0f af d0          	imul   %r8,%rdx
   1bc03:	48 29 d7             	sub    %rdx,%rdi
   1bc06:	4a 8d 14 07          	lea    (%rdi,%r8,1),%rdx
   1bc0a:	48 39 f8             	cmp    %rdi,%rax
   1bc0d:	48 0f 46 fa          	cmovbe %rdx,%rdi
   1bc11:	48 89 f8             	mov    %rdi,%rax
   1bc14:	4c 29 c0             	sub    %r8,%rax
   1bc17:	4c 39 c7             	cmp    %r8,%rdi
   1bc1a:	48 0f 42 c7          	cmovb  %rdi,%rax
   1bc1e:	c3                   	retq   
*/
