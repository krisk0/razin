RAZIN - Rigorous algebraic zero impeding numeration
=====

For now, I offer
  a) A subroutine to calculate HNF of a matrice whose determinant is small in absolute value;
  b) PLU decomposition and inversion of a matrice modulo 2^64;
  c) Inversion of an integer upper-triangular matrice with small determinant and positive diagonal;
  d) Python binding to some methods of FLINT, called flint_sage;
  e) A subroutine to count an integer matrice determinant.
 
Where are the files?
^^^^^^^^^^^^^^^^^^^^
See *filelist.txt*

My new algorithm to compute integer matrice determinant: ``fmpz_mat_det_hermitian_decomposition()``

Small-det HNF: ``mod_mat_HNF_nonsquare()`` and ``mod_mat_HNF()``

What's the target?
^^^^^^^^^^^^^^^^^^
Mid-range target of the project is fast HNF computation with a new algorithm inspired by W.Stein double-determinant. Python wrapper is minimalistic and only contains functions required to reach the goal or to test/benchmark subroutines/algorithms.

The fastest in the open-source world subroutines to count big matrice determinant and small-det matrice Hermite form are a by-product.

FLINT and GMP data structures and subroutines are used.

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
