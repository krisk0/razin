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
            
            test if fmpz_mat_fflu, fmpz_mat_rref, fmpz_mat_rref_fraction_free
             do smth useful for me. If they count HNF, benchmark against NTL
             HNF
