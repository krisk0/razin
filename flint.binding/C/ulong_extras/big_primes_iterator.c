// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include "ulong_extras_.h"
#include "is_big_sureprime.c"

#if USE_BASE7_MR
 #define MIN_n_primes_rev 299210839
#else
 #define MIN_n_primes_rev 103
#endif

// Code in this file in known to work on amd64. Don't know what will happen
//  on other arch

// Structure used to quickly iterate thru primes in descending order, once or
//  many times. Max allowed range: MIN_n_primes_rev..2**64-59 where 2**64-59 is 
//  maximal native prime and MIN_n_primes_rev=299210839 is minimal number for
//  which trial_gcd_test_then_big_sureprime_test() works correctly

#define LOUD_ALLOC 1

#if (defined (__amd64__) || defined (__i386__) || defined (__i486__)) 
 #define GCD n_gcd_odd_odd
 mp_limb_t n_gcd_odd_odd(mp_limb_t x,mp_limb_t y);
#else
 #define GCD n_gcd_full
#endif
#define GCDt(x) if( 1<GCD( n, WORD(x) )) return 0;

int
trial_gcd_test_then_big_sureprime_test(mp_limb_t n)
 {
  GCDt( 0xF1354E62564E313  ) //7*11*13*17*19*23*29*31*37*41*43*47*53
  GCDt( 0x6329899EA9F2714B ) //59*61*67*71*73*79*83*89*97*101
  GCDt( 0x21A3907D1B750A13 ) //193*407521*299210837*103
  // adding line below did not give any noticeable time difference
  GCDt( 0x825F18A4856CE7FB ) //107*109*113*127*131*137*139*149*151
  return n_is_big_sureprime(n);
 }

int trial_gcd_test_then_BPSW(mp_limb_t n)
/*
Presumably n_is_big_sureprime() and n_is_probabprime_BPSW() never make 
 mistakes 
*/
 {
  GCDt( 0xF1354E62564E313  ) //7*11*13*17*19*23*29*31*37*41*43*47*53
  GCDt( 0x6329899EA9F2714B ) //59*61*67*71*73*79*83*89*97*101
  // adding GCDt below slows down a little, so disabling it
  //GCDt( 0x58EDCB4C9ED39C8B ) //103*107*109*113*127*131*137*139*149
  return n_is_probabprime_BPSW(n);
 }

#undef GCDt
#undef GCD

int prevmod30[]={ 1, 2, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 1, 2, 1, 2,
                  3, 4, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6 };

static __inline__ void
n_primes_rev_shift(n_primes_rev_t i)
// shift state and get ready to output next number
 {
  // Are we done?
  mp_limb_t curr=i->numbers[i->index];
  if(curr == 1)
   return;
  // Is this not first time?
  slong curr_mod_30=i->last_output_mod_30;
  if(curr_mod_30 < 0)
   {
    if(i->index+1 < i->allocated_size )
     {
      // re-use already found number
      ++i->index;
      return;
     }
    curr_mod_30 = curr % 30;
   }
  // Is this even number?
  int grow_index=curr & 1;
  // set curr to previous prime, or 1 if range exhausted
  while(1)
   {
    // curr==prime already output or even number, curr_mod_30==curr % 30
    mp_limb_t shi=prevmod30[ curr_mod_30 ];
    curr -= shi;
    if(curr < MIN_n_primes_rev)
     {
      curr = 1;
      break;
     }
    curr_mod_30 -= (slong)shi;
    curr_mod_30 += 30*(curr_mod_30<0);
    #if USE_BASE7_MR
     if( trial_gcd_test_then_big_sureprime_test(curr) )
    #else
     if( trial_gcd_test_then_BPSW(curr) )
    #endif
      break;
   }
  if(grow_index)
   {
    // if previous number was odd, advance index
    if( ++i->index >= i->allocated_size )
     {
      #if LOUD_ALLOC
       gmp_printf("n_primes_rev_shift(): %d==%d, old ptr=%MX",i->index,
        i->allocated_size,(mp_limb_t)i->numbers);
      #endif
      i->numbers=realloc( i->numbers, 
       sizeof(mp_limb_t)*(i->allocated_size += 10) );
      #if LOUD_ALLOC
       gmp_printf(", reallocated %d, new ptr=%MX\n",i->allocated_size,
        (mp_limb_t)i->numbers);
      #endif
     }
    #if LOUD_ALLOC
     else
      {
       flint_printf("n_primes_rev_shift(): %d<%d\n",i->index,i->allocated_size);
      }
    #endif
    i->numbers[ i->index ] = curr;
   }
  else
   // if previous number was even, replace it
   i->numbers[ i->index ] = curr;
  i->last_output_mod_30=curr_mod_30;
 }

mp_limb_t 
n_primes_rev_reset(n_primes_rev_t i)
// un-iterate back, to re-iterate from the start
 {
  #if LOUD_ALLOC
   flint_printf("n_primes_rev_reset() setting allocated_size to %d\n",i->index+1);
  #endif
  i->allocated_size=i->index+1;
  i->index=0;
  i->last_output_mod_30=-33;
  return i->numbers[0];
 }

// TODO: make it possible to pre-allocate array of custom length rather than 10
mp_limb_t
n_primes_rev_init(n_primes_rev_t i,mp_limb_t stArt)
/*
initialize iterator i
 
stArt is odd prime > MIN_n_primes_rev
      or 0 (which means start from maximal 64-bit prime 2**64-59)
      or even number (which means start from maximal prime < stArt)

return the first prime
*/
 {
  i->index=0;
  i->allocated_size=10;
  i->numbers=flint_malloc( 10 * sizeof(mp_limb_t) );
  #if LOUD_ALLOC
   gmp_printf("allocated_size initialized with 10, ptr=%MX\n",(mp_limb_t)i->numbers);
  #endif
  if( !stArt )
   stArt=WORD( 0xFFFFFFFFFFFFFFC5 );
  i->numbers[0] = stArt;
  i->last_output_mod_30=stArt % 30;
  if( stArt & 1 )
   return stArt;
  n_primes_rev_shift(i);
  return i->numbers[0];
 }

void 
n_primes_rev_clear(n_primes_rev_t i)
 {
  #if LOUD_ALLOC
   flint_printf("n_primes_rev_clear()\n");
  #endif
  flint_free(i->numbers);
 }

mp_limb_t
n_primes_rev_show_again(n_primes_rev_t i)
/*
re-output last number. Useful if n_primes_rev_init() or n_primes_rev_reset()
 was called in a different subroutine
*/
 {
  return i->numbers[i->index];
 }

mp_limb_t
n_primes_rev_next(n_primes_rev_t i)
/*
return previous prime down to MIN_n_primes_rev, 
 or 1 if range exhausted
*/
 {
  n_primes_rev_shift(i);
  return i->numbers[i->index];
 }

mp_limb_t
count_primes_in_range(mp_limb_t lo,mp_limb_t up)
/*
return count of primes in range lo..up, up=0 means FFFFFFFFFFFFFFC5
lo should be in range MIN_n_primes_rev..2**64
*/
 {
  n_primes_rev_t i;
  mp_limb_t x=n_primes_rev_init(i,up);
  mp_limb_t co=0;
  while(1)
   {
    if(x<lo)
     {
      n_primes_rev_clear(i);
      return co;
     }
    ++co;
    x=n_primes_rev_next(i);
   }
 }

void
loud_primes_count(mp_limb_t lo)
 {
  gmp_printf("count of primes not smaller than 0x%MX\n",lo);
  mp_limb_t c=count_primes_in_range(lo,0);
  gmp_printf("... = 0x%MX\n",c);
 }

#undef LOUD_ALLOC
