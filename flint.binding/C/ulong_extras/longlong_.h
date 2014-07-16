#ifndef LONGLONG__H
#define LONGLONG__H

#if (GMP_LIMB_BITS == 64 && defined (__amd64__))

#define add_sssaaa0aa(sh, sm, sl, ah, am, al,    bm, bl)  \
  __asm__ ("addq %8,%q2\n\tadcq %6,%q1\n\tadcq %4,%q0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "N"   (0x0),              \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  \

#endif

#endif
