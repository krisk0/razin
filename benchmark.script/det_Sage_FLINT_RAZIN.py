#!/usr/bin/python2 -u
# -*- coding: utf-8

# This program is part of RAZIN
# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# This program requires flint_sage package

'''
Benchmark determinant calculation: three Sage subroutines, FLINT, RAZIN block 
 Gauss
'''

import sage.all,numpy,time,sys,os,signal
import flint_sage as flint
'''
standard ways to limit execution time do not work on Sage .determinant(), so use
 Sage .Fork(timeout=...)
'''
import sage.parallel

dim=300
experiments=3
ways_desc='Sage NTL pari FLINT 4blck'
ways=5
entry_max=1<<1000
time_constraint=time_limit=None        # to be detected auto-magicallly
loud=os.environ.get('loud')
sage.all.set_random_seed('20140714')
ZZ=sage.all.ZZ

t1=numpy.resize( numpy.array( [], dtype=float ), ways )

def signal_handler(signum, frame):
 '''
 stackoverflow.com/questions/366682/
  how-to-limit-execution-time-of-a-function-call-in-python
 '''
 pass

signal.signal(signal.SIGALRM, signal_handler)

def random_matrice():
 return sage.all.random_matrix(ZZ, dim, x=-entry_max, y = entry_max)

def count_det(a,me,i):
 global t1,time_limit,time_constraint
 if t1[i]<-0.1:
  return
 if loud:
  print 'starts method %s' % me
 if time_limit != None:
  F=time_limit(a.determinant)
  t0=time.time()
  try:
   one=F(algorithm=me,proof=True)
   t1[i] += time.time()-t0
   if loud:
    print 'spent',t1[i]
  except:
   # why "[Errno 4] Interrupted system call" is printed? I did not ask that
   t1[i]=-1
   if loud:
    print 'TiMeOuT'
 else:
  '''
   time_constraint magic: set it 30% more than time spent by Sage padic 
    method
  '''
  t0=time.time()
  one=a.determinant(algorithm=me,proof=True)
  t0=time.time()-t0
  t1[i] += t0
  if loud:
   print 'spent',t0
  time_constraint=int( t0*1.3 ) # signal.alarm wants integer?
  time_limit=sage.parallel.decorate.Fork(timeout=time_constraint)

def det_4_ways(a):
 way=0
 # 'linbox' is the 4th way, disqualified for providing no sure way
 for m in 'padic','ntl','pari': 
  count_det(a,m,way)
  a._clear_cache()
  way += 1

def det_flint(a,f,way):
 global t1
 if t1[way]<-0.1:
  return
 aF=flint.fmpz_mat( a )
 if loud:
  print 'starts subroutine %s' % f
 t0=time.time()
 signal.alarm(time_constraint)
 try:
  one=f(aF)
  t1[way] += time.time()-t0
  if loud:
   print 'spent',t1[way]
 except:
  t1[way]=-1
  if loud:
   print 'TiMeOuT'

def chew_desc( d ):
 r,l=d.split(' '),0
 for x in r:
  s=len(x)
  if s>l:
   l=s
 return r,l

def pretty_print_result():
 print
 desc,max_l=chew_desc( ways_desc )
 print (' '*max_l)+'    time'
 print ('-'*max_l)+' ---------'
 fmt='%'+str(max_l)+'s  %.2e'
 fm2='%'+str(max_l)+'s  Time-out'
 for i in range(ways):
  if t1[i]>0:
   print fmt % (desc[i],t1[i])
  else:
   print fm2 % desc[i]

for x in range(experiments):
 a=random_matrice()
 det_4_ways(a)
 det_flint(a,flint.det,3)
 det_flint(a,flint.det_20140704,4)
 for x in range(ways):
  sys.stdout.write( '%.4e ' % t1[x] )
 print

pretty_print_result()
'''
with stupid naive vector dot multiplication results were:
 experiments=3 dim=200:
          time
 ----- ---------
  Sage  2.00e+02
   NTL  1.39e+02
  pari  Time-out
 FLINT  7.38e+01
 4blck  6.82e+01
 
 experiments=3 dim=250:
          time
 ----- ---------
  Sage  3.06e+02
   NTL  3.11e+02
  pari  Time-out
 FLINT  1.60e+02
 4blck  1.22e+02
 
 experiments=3 dim=300:
          time
 ----- ---------
  Sage  4.34e+02
   NTL  Time-out
  pari  Time-out
 FLINT  3.04e+02
 4blck  1.99e+02
'''
