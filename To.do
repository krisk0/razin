11 Mar 2014 Add Python interface to FLINT functions solving matrice 
             equation and counting HNF over Z/nZ.
             
            fmpz_mat_solve+ nmod_mat_rref 
             
            Test the new functions by implementing a fragment of Stein
             double-determinant algorithm. When test passes, benchmark.

12 Mar 2014 We don't want slow procedure when faster is available
            
            Before proceeding to nmod_mat_rref, implement a wrapper to 
             fmpq_mat_solve_dixon(). And benchmark

            fmpq_mat constructor+    fmpq_mat_solve_dixon+ fmpq_mat_inv+
            
            benchmark Dixon/solve_right with one vector
            
            convert fmpq_mat to Sage equivalent *not done* Do I need it?
