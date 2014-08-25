11 Mar 2014 Add Python interface to FLINT functions solving matrice 
             equation and counting HNF over Z/nZ.
             
            fmpz_mat_solve+ nmod_mat_rref --- ups, this procedure does not
             work for non-prime moduli
             
            Test the new functions by implementing a fragment of Stein
             double-determinant algorithm. When test passes, benchmark.

12 Mar 2014 We don't want slow procedure when faster is available
            
            Before proceeding to nmod_mat_rref, implement a wrapper to 
             fmpq_mat_solve_dixon(). And benchmark --- done

            fmpq_mat constructor+    fmpq_mat_solve_dixon+ fmpq_mat_inv+
            
            benchmark Dixon/solve_right with one vector --- done
            
            convert fmpq_mat to Sage equivalent *not done* Do I need it?
            
15 Mar 2014 test if fmpz_mat_fflu, fmpz_mat_rref, fmpz_mat_rref_fraction_free
             do smth useful for me. If they count HNF, benchmark against NTL
             HNF. 

            --- Done. The procedures return strange result, destroy input
            matrice, and fmpz_mat_rref_fraction_free is no longer there.

            Bill Hart confirmed that FLINT currently does not do row 
            transformation on matrice over residue ring. So I will implement
            nmod_mat_HNF() for square non-singular matrice.
             
            
17 Mar 2014 nmod_mat_HNF() test passes, now benchmark. What subroutines to 
             benchmark against except Sage _hnf_modn() and NTL mat_ZZ.h HNF?

            benchmark_nmod_mat_HNF-sage.py says my beautiful C code is slower
             than ugly Cython ._hnf_modn_impl()
             
            where did I make mistake? Arrr, I am doing 4 scalar multiplications
             even if pivot element is 1
            
            Rewrite nmod_mat_HNF-debug.c. Introduce 3 subroutines instead of one
             (one for general case, one for small modulo, one for degree of 2)
             --- delayed

19 Mar 2014 nmod_mat_HNF-debug.c is at least 2 times faster than Sage method.
            // 402 lines/4 days = 100 lines of code per day. Awfully slow

            To further speed up, fix subroutine DKryskov_nmod_early_abort() and
             specialize the code for small modulo (<2**32). Also can shorten 
             some vectors by 1 when doing row operations
             
            The task in paragraph above delayed for indefinite period

20 Mar 2014 re-implement Stein double-determinant algorithm (proof=1) with no 
             modification, just using faster subroutines for salvation of linear
             system and calculating HNF of matrice with small determinant --- done
            
            profile: see what fragment takes most times

22 Mar 2014 found Strassen multiplication algorithm in Matvei Nazaruk Git

            https://github.com/matveinazaruk/Strassen/commit/
             ef4eda1de8ce193a47204fa67ec42862f579b54e

            if the code compiles, maybe include it in Razin? --- delayed

            After finding H = HNF(A) I need to find transformation matrice U such
             that U * H = A.
            Obviously no FLINT function exists that does it efficiently. Should I
             write my own?
            Ups, I can find that U, but it won't do any good because HNF computation
            is very fast compared to other matrice operation. For instance
            computing det U is terribly slow. Wrong way, go back.

24 Mar 2014 Profile the double-det algorithm: find out time spent by
             Dixon linear solver,
             solve_system_with_difficult_last_row(),
             add_col(),
             1st add_row(),
             2nd add_row(),
             double_det(),
             det_given_divisor() --- profiling done

            record g values and h.b.(A)/h.b.(U) where A and U
            are as above and h.b.() is Hadamard bound on determinant;
            record count of tries in solve_system_with_difficult_last_row() 
             --- delayed, got better idea

26 Mar 2014 arrr, Dixon / IML efficiency in solve_system_with_difficult_last_row()
             floats from test to test. Must do more tests, this takes time, 
             especially for n=1000

27 Mar 2014 To speed up solve_system_with_difficult_last_row(), need a 
             subroutine to quickly put unimodular matrice to row-echelon form
             and find transformation that does it

03 Apr 2014 tmod_mat_PLU_mod_machine_word() too slow, need to write recursive
             version
             
            I wonder why time of solve_system_with_difficult_last_row() does 
             not equal sum of times: for n=3000, bits=8 
             
            752.7781+792.1059 = 1544.8 < 1671.9109     that's 8% losses 

            Will investigate this

04 Apr 2014 Need to
             1. Improve check of tmod_mat_solver() in test_tmod_mat.py
              --- done
             2. Implement recursive algorithm for inversion of lower-tria matrix 
              --- done
             3. Benchmark the algorithm including LU decomposition and inverse 
                 for big and fat matrices, see how big diagonal entries are
              --- done, tmod_mat is fast enough, so I delay implementing 
              recursive version of tmod_mat_PLU_mod_machine_word()

06 Apr 2014 Test PLU decomposition for square input ??? do I need it?

11 Apr 2014 After MPIR-2.7.0 test terminates, install GMP properly and re-test.

            Finally, write a code to replace Stein double-det.

25 Aug 2014 ? or better speed-up my fastest in the open-world determinant 
             calculation?
