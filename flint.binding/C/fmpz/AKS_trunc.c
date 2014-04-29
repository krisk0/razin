// This program is part of RAZIN
// Copyright 袁轶君 (Yijun Yuan), Денис Крыськов 2014
// Licence: public domain

// original code was not working for large n, fixes by Денис Крыськов make it work
// the fixed lines are marked `Денис Крыськов was here`

// This file contains public domain code put online by 袁轶君 (Yijun Yuan) 23 apr
//  2014 at https://github.com/YijunYuan/AKS_FLINT

#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fmpz_mod_poly.h>
#include<flint/arith.h>
//#include<mpir.h>
#include<stdio.h>
#include<math.h>

static __inline__ void 
MULTIPLICATIVE_ORDER(fmpz_t out,fmpz_t n,fmpz_t k){
    fmpz_gcd(out,n,k);
    if(!fmpz_equal_ui(out,1)){
        fmpz_set_si(out,-1);return;
    }
    fmpz_t i;fmpz_init(i);
    for(fmpz_one(i);;fmpz_add_ui(i,i,1)){
        fmpz_powm(out,k,i,n);
        if(fmpz_equal_ui(out,1)){
                fmpz_set(out,i);
                fmpz_clear(i);
                return;
        }
    }
}

int 
AKS(fmpz_t n)
// n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64
{
 // Денис Крыськов was here - some lines deleted in view of specification
    fmpz_t temp1;fmpz_init(temp1);
 /*Step1*/
    mpz_t k;mpz_init(k);fmpz_get_mpz(k,n);
 if(mpz_perfect_power_p(k)!=0){
            fmpz_clear(temp1);mpz_clear(k);
            return 0;
    }
 mpz_clear(k);
 /*Step2*/
 fmpz_t c,r;fmpz_init(c);fmpz_init(r);
 fmpz_set_ui(c,floor(sqrt(fmpz_dlog(n))));//c=[sqrt(log(n))]
    for(fmpz_set_ui(r,2);;fmpz_add_ui(r,r,1)){
        MULTIPLICATIVE_ORDER(temp1,r,n);
        if(fmpz_cmp(temp1,c)>0)break;
    }
    /*Step3*/
    fmpz_t a;fmpz_init(a);
    for(fmpz_one(a);fmpz_cmp(a,r)<=0;fmpz_add_ui(a,a,1)){
        fmpz_gcd(temp1,a,n);
        if(fmpz_cmp(temp1,n)<0&&fmpz_cmp_ui(temp1,1)>0){
            fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
            return 0;
        }
    }
    /*Step4*/
    if(fmpz_cmp(n,r)<=0){
            fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
            return 1;
    }
    /*Step5*/
    fmpz_set_d(c,fmpz_dlog(n));//c=log(n)
    arith_euler_phi(temp1,r);//temp1=euler(r)
    fmpz_mul(c,c,temp1);//c=log(n)*euler(r)
    fmpz_clear(temp1);
    fmpz_mod_poly_t e;fmpz_mod_poly_init(e,n);
    mp_limb_t n_ui=fmpz_get_ui(n);
    mp_limb_t r_ui=fmpz_get_ui(r);
    for(fmpz_one(a);fmpz_cmp(a,c)<=0;fmpz_add_ui(a,a,1)){
        fmpz_mod_poly_t modulo;fmpz_mod_poly_init(modulo,n);
        fmpz_mod_poly_t p     ;fmpz_mod_poly_init(p     ,n);
        fmpz_mod_poly_t q     ;fmpz_mod_poly_init(q     ,n);
        fmpz_mod_poly_set_coeff_ui(modulo,0             ,n_ui-1); // Денис Крыськов was here
        fmpz_mod_poly_set_coeff_ui(modulo,r_ui          , 1);//modulo=x^r-1
        fmpz_mod_poly_set_coeff_ui(p     ,1             , 1);
        fmpz_mod_poly_set_coeff_fmpz(p   ,0             , a);//p=x+a
        //fmpz_mod_poly_set_coeff_ui(q     ,fmpz_get_ui(n), 1);
        // The commented line above is not good for big n, if below fixes it
        if(n_ui > r_ui) // Денис Крыськов was here
         fmpz_mod_poly_set_coeff_ui(q     ,n_ui % r_ui, 1);
        else
         fmpz_mod_poly_set_coeff_ui(q     ,n_ui       , 1);
        fmpz_mod_poly_set_coeff_fmpz(q   ,0             , a);//q=x^n+a or x^(n-r)+a
        fmpz_mod_poly_powmod_fmpz_binexp(p,p,n,modulo);//p=p^n mod modulo
        fmpz_mod_poly_divrem(e,q,q,modulo);//q=q mod modulo
        if(fmpz_mod_poly_equal(p,q)==0){
                fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
                fmpz_mod_poly_clear(modulo);fmpz_mod_poly_clear(p);fmpz_mod_poly_clear(q);fmpz_mod_poly_clear(e);
                return 0;
        }
        fmpz_mod_poly_clear(modulo);fmpz_mod_poly_clear(p);fmpz_mod_poly_clear(q);
    }
    fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
    fmpz_mod_poly_clear(e);
 return 1;
}

int 
AKS_ui(mp_limb_t n)
// n must be odd in range 3..ULONG_MAX where ULONG_MAX equals 2**64-1 for amd64
 {
  fmpz_t m; fmpz_init_set_ui( m, n );
  int rc=AKS(m);
  fmpz_clear(m);
  return rc;
 }
