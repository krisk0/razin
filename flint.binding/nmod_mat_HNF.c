#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

slong
nmod_mat_HNF(nmod_mat_t A)
/*
HNF: Hermite Normal Form as defined by W.Stein/C.Pernet

A: non-singular over integers, det(A) is a divisor of n=mod
let H be HNF of A over Z/nZ

if H is non-singular, modify A to be equal H and return -1
if H is singular, return j such that H[j,j] is zero, and return A with
 the property

 0 <= A[i,j] < n,  A[i,j] = H[i,j] modulo n for i<j
 (all other entries except j-1 defined by line above undefined)
*/
 {
 }
