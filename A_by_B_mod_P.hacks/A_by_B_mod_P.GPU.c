// Put into public domain by Sardar Anisul Haque, Marc Moreno Maza
// University of Western Ontario, London N6A 1M8, Canada
// works for p ?not larger 30 bits?

__device__ int double_mul_mod(int a, int b, int p, double pinv) 
 {
  int q = (int) ((((double) a) * ((double) b)) * pinv);
  int res = a * b - q * p;
  return (res < 0) ? (-res) : res;
 }
 
