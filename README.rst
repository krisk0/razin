RAZIN - Rigorous algebraic zero impeding numeration
=====

For now, I offer
 
  a) A subroutine to calculate HNF of a square integer matrice with low determinant;
  b) Python binding to some methods of FLINT;
  c) some algorithms benchmarks (such as determinant, linear equations solver)
 
Where are the files?
^^^^^^^^^^^^^^^^^^^^
See *filelist.txt*. The fastest in galaxy subroutine to compute HNF via modular technique is in flint.binding/nmod_mat_HNF.c

What's the target?
^^^^^^^^^^^^^^^^^^
Mid-range target is fast HNF computation with an algorithm resembling W.Stein double-determinant. Python wrapper is minimalistic and only contains functions required to reach the goal or to test/benchmark subroutines/algorithms.

Only FLINT functions will be used for solving hard sub-problems, unless a function from another library (NTL, IML, LinBox, ...) turn to be faster or a serious problem with FLINT discovered. Presently components of FLINT wrapped into flint_sage Python package work as they should. If you think otherwise, your bug-report is welcome.

*I wrote the paragraph above in the morning, to discover in the evening that FLINT is not so good at inverting big matrices as dear old Sage*

Glossary
^^^^^^^^

:FLINT:
    C numerical/matrice library with save-every-penny approch to arithmetic. `Here <http://www.flintlib.org/>`_

:krisk0/razin:
    a project containing Python/Cython code for integer linear algebra attempting to use the fastest available algorithms and save every penny

:NTL:
    C++ numerical/matrice library with gluttonous bigint constructor. `Here <http://shoup.net/ntl/>`_

:Captain Flint: 
    fictional Caribbean adventurer (who said pirate?)
nmod_mat_HNF.c
:Stepan Razin: 
    16??-1671, kozak, rebel leader (who said gangster?)
