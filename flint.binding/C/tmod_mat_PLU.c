// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

long 
tmod_mat_PLU_mod_machine_word(long* PR,tmod_mat_t S)
/*
attempt to find PLU factorisation of S modulo 2**64 on amd64, such that
 P*original S = L * U

lower line of permutation P stored into PR, followed by array R such that
R[i] * diagonal(U)[i] = 1 modulo 2**64

return 1 on success, 0 on failure

S modified. On success S contains LU in compressed form, just like FLINT 
*/
 {
 }
