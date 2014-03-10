RAZIN - Rigorous algebraic zero impeding numeration
=====

For now, nothing big except
 
 a) Python binding to FLINT integer matrice;
 b) Unimodularity test benchmarks::
 
                              U      warped U
                     Sage  2.13e+02  2.00e+02
                      NTL  6.13e-05  3.05e+01
                     pari  5.32e-05  4.50e+01
                    FLINT  2.29e+01  2.30e+01
                  solve_r  1.23e+03  7.39e+00

Benchmark explained: 
1) if matrice is unimodular, use PARI/GP or NTL (Sage native subroutine 4 million times slower).
2) if matrice is warped unimodular, FLINT is slightly faster (32% faster than NTL).

Yes, Python is slower than C, but not 4e6 times! Could it be that NTL and PARI/GP determinant calculation is broken? Any example?

My e-mail address is at the bottom of "Hello world" page of my `blog <http://tiny.cc/DKryskov>`_.

Where are the files?
^^^^^^^^^^^^^^^^^^^^
See *filelist.txt*

Todo
^^^^
Quick integer linear algebra subroutines based on FLINT will be here s0oner or later. 

Glossary
^^^^^^^^
:FLINT:
    C numerical/matrice library with save-every-penny approch to arithmetic. `Here <http://www.flintlib.org/>`_
:NTL:
    C++ numerical/matrice library with gluttonous bigint constructor but clever algorithms. `Here
    <http://shoup.net/ntl/>`_
:Captain J. Flint: 
    fictional Caribbean adventurer (who said pirate?)
:Stepan Razin: 
    16xx-1671, kozak, rebel leader (who said bandit?)
