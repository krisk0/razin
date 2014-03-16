#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

slong
nmod_mat_HNF(nmod_mat_t A)
/*
HNF: Hermite Normal Form as defined by W.Stein/C.Pernet

Input data:
~~~~~~~~~~
A: square non-singular over integers, A.r>0, det(A) is a divisor of n=mod

let H be HNF of A over Z/nZ

Output data:
~~~~~~~~~~~
if H is non-singular, modify A to be equal H and return -1
if H is singular, return j such that H[j,j] is zero, and return A with
 the property
 for all i in 0..j-1 ( 0 <= A[i,j] < n,  A[i,j] = H[i,j] modulo n )
 (all other entries of A except j-1 defined by line above undefined)
*/
 {
 }
