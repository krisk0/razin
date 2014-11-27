// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#include <stdio.h>
#include <string.h>
#include <flint/flint.h>

char* 
GMP_or_MPIR_version_c()
/*
GMP or MPIR lib called libgmp* or libmpir* should define gmp_version and/or 
 mpir_version
 
this subroutine extracts this string (using malloc to allocate memory)
*/
 {
  slong s;
  char* q;
  const char* r;
  #if defined(__MPIR_VERSION)
   s=strlen(mpir_version)+5;   
   r=mpir_version;
   #define p "MPIR"
  #else
   s=strlen(gmp_version)+4;
   r=gmp_version;
   #define p "GMP"
  #endif
  q=malloc(s);
  sprintf(q,"%s %s",p,r);
  return q;
 }
#undef p
