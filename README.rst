RAZIN - Rigorous algebraic zero impeding numeration
=====

For now, I only offer
 
 a) Python binding to FLINT integer matrice;
 b) Unimodularity test benchmarks. **FLINT fastest, NTL on 2nd place**
 c) Solve_right() benchmarks. **fmpz_mat_solve() tested and benchmarked**, mpq_mat_solve_dixon() will be tested soon. The latter is expected to be faster on big matrix than the former.
 
Where are the files?
^^^^^^^^^^^^^^^^^^^^
See *filelist.txt*. 

What's the target?
^^^^^^^^^^^^^^^^^^
Mid-range target is fast HNF computation with an algorithm resembling W.Stein double-determinant. Python wrapper is minimalistic and only contains functions required to reach the goal or to test/benchmark/experiment with algorithms.

Only FLINT functions will be used for solving hard sub-problems, unless functions from another library (NTL, IML, LinBox, ...) turn to be faster or a serious problem with FLINT discovered. Presently components of FLINT wrapped into flint_sage Python package work as they should. If you think otherwise, your bug-report is welcome.

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

:Stepan Razin: 
    16??-1671, kozak, rebel leader (who said gangster?)
