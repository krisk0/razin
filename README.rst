RAZIN - Rigorous algebraic zero impeding numeration
=====

For now, nothing big except
 
 a) Python binding to FLINT integer matrice;
 b) Unimodularity test benchmarks. 
 
5 algorithms testing big matrice for unimodularity were benchmarked, results in the table below::
 
                              U      warped U
                  ------- --------- ---------
                     Sage  2.14e+02  1.98e+02
                      NTL  3.01e+01  2.98e+01
                     pari  4.51e+01  4.01e+01
                    FLINT  2.26e+01  2.27e+01
                  solve_r  1.23e+03  7.34e+00

Benchmark result verdict: 

1) if matrice is unimodular, use FLINT (NTL on 2nd place, 33% slower).
2) if matrice is warped unimodular, solve_right() quickly discovers it.

FLINT counts determinant faster than others!

Where are the files?
^^^^^^^^^^^^^^^^^^^^
See *filelist.txt*. How data is generated and what is benchmarked? See *unimodular.5way.cout* and *unimodular.5way.py*

Todo
^^^^
Quick integer linear algebra subroutines based on FLINT will be here s0oner or later. 

Glossary
^^^^^^^^

:FLINT:
    C numerical/matrice library with save-every-penny approch to arithmetic. `Here <http://www.flintlib.org/>`_

:NTL:
    C++ numerical/matrice library with gluttonous bigint constructor. `Here <http://shoup.net/ntl/>`_

:Captain Flint: 
    fictional Caribbean adventurer (who said pirate?)

:Stepan Razin: 
    16??-1671, kozak, rebel leader (who said gangster?)
