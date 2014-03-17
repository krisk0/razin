11 Mar 2014 Add Python interface to FLINT functions solving matrice 
             equation and counting HNF over Z/nZ.
             
            fmpz_mat_solve+ nmod_mat_rref --- ups, this procedure does not
             for non-prime moduli
             
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
             
16 Mar 2014 --- Done. The procedures return strange result, destroy input
             matrice, and fmpz_mat_rref_fraction_free is no longer there.
             
            Bill Hart confirmed that FLINT currently does not do row 
            transformation on matrice over residue ring. So I will implement
            nmod_mat_HNF() for square non-singular matrice.
            
17 Mar 2014 Test passes, now benchmark. What subroutines to benchmark against
             except Sage _hnf_modn() and NTL mat_ZZ.h HNF?
