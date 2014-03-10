#include <iostream>

__uint64_t mulmod1(__uint64_t a, __uint64_t b, __uint64_t m) {
  if (b < a)
    std::swap(a, b);
  __uint64_t res = 0;
  for (__uint64_t i = 0; i < a; i++) {
    res += b;
    res %= m;
  }
  return res;
}

__uint64_t mulmod2(__uint64_t a, __uint64_t b, __uint64_t m) {
  __uint64_t r;
  __asm__
  ( "mulq %2\n\t"
      "divq %3"
      : "=&d" (r), "+%a" (a)
      : "rm" (b), "rm" (m)
      : "cc"
  );
  return r;
}  
#if 0
 cc clobber conditional code register
 =& early clobber - dx cannot be input register
 + input and output
 % can be switched with next operand
#endif

int main() {
  using namespace std;
  __uint64_t a = 853467ULL;
  __uint64_t b = 21660421200929ULL;
  __uint64_t c = 100000000000007ULL;

  cout << mulmod1(a, b, c) << endl;
  cout << mulmod2(a, b, c) << endl;
  return 0;
}
