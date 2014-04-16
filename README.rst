RAZIN - Rigorous algebraic zero impeding numeration
=====

For now, I offer
  a) Two subroutines to quickly calculate HNF of an integer matrice in two practically important cases:
         a) it is non-singular;
         b) its determinant is known and lies in range 2..2^64-1 where 2^64-1 is maximal native unsigned integer;
  b) PLU decomposition of a matrice modulo 2^64;
  c) inversion of square lower-triangular matrice modulo 2^64;
  d) Python binding to some methods of FLINT, called flint_sage;
  e) some algorithms benchmarks (such as determinant, linear equations solver)
 
Where are the files?
^^^^^^^^^^^^^^^^^^^^
See *filelist.txt*.

The fastest in open-source world subroutines to compute HNF of a matrice with small in absolute value determinant are callable from C or Python; C subroutines are called ```mod_mat_HNF_nonsquare()``` and ```mod_mat_HNF()```; see file *flint.binding/nmod_mat_HNF.c*

W.Stein double-deteminant algorithm to compute HNF specialized for non-singular matrice and instructed to sometimes use faster low-level routines like FLINT Dixon lifting instead of Sage method is in *benchmark.script/profile_Sage_hnf_square.py*. The modified procedure is faster than original, see bottom of *profile-dd_algorithm.cout* for details.

PLU decomposition modulo 2^64 test: *flint.binding/test_tmod_mat.py*

What's the target?
^^^^^^^^^^^^^^^^^^
Mid-range target of the project is fast HNF computation with a new algorithm inspired by W.Stein double-determinant. Python wrapper is minimalistic and only contains functions required to reach the goal or to test/benchmark subroutines/algorithms.

Only FLINT functions will be used for solving hard sub-problems, unless a function from another library (NTL, IML, LinBox, ...) turn to be faster or a serious problem with FLINT discovered. Presently components of FLINT wrapped into flint_sage Python package work as they should. If you think otherwise, your bug-report is welcome.

*I wrote the paragraph above in the morning, to discover in the evening that FLINT is not so good at inverting big matrices as dear old Sage*

Bug?
^^^^
If you find out that some subroutine of Razin produces bad result or crashes, your bug-report is welcome. If it includes the following data: sample code including input data, description of what happens, expected result (in case you think output is wrong). For instance, if you think that fmpz_mat_hermite_form() works incorrectly, provide
  a) your code forming input matrice,
  b) output matrice,
  c) the correct output

Glossary
^^^^^^^^

:FLINT:
    C numerical/matrice library with save-every-penny approach to arithmetic. `Here <http://www.flintlib.org/>`_

:krisk0/razin:
    a research log and software package with quickest subroutines for integer linear algebra 

:NTL:
    C++ numerical/matrice library with gluttonous bigint constructor. `Here <http://shoup.net/ntl/>`_

:Captain Flint: 
    fictional Caribbean adventurer (who said pirate?)

:Stepan Razin: 
    16??-1671, kozak, rebel leader (who said gangster?)
