# -*- coding: utf-8

# This program is part of RAZIN

# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# This program locates sage library, then builds flint_sage Python module

# If this program fails to find sage include, set MY_SAGE_IS_HERE to point to 
#  directory containing sage/rings/integer

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os,sys
osE=os.environ.get
write=sys.stdout.write

def find_sage_include_dir( prefix, suffix, user_set ):
 if user_set == None:
  # no user-supplied directory, start with empty list
  w=[]
 else:
  # user-supplied setting tried on the 1st place
  w=[user_set]
 if prefix==None:
  # EPREFIX unset, only one variant
  w.append( suffix )
 else:
  # EPREFIX set, try both ways
  w += [ prefix+suffix, suffix ]
 r=None
 for x in w:
  if os.path.isdir(x):
   r=x
   break
 if r==None:
  # todo: use locate command
  write('Failed to find include dir, tried directories:\n')
  for x in w:
   write('   >'+x+'<\n')
  sys.exit(1)
 write('sage lives in '+r+'\n')
 return r

'''
 sage-on-gentoo ebuilds put *.p* into $EPREFIX/usr/share/sage
'''

p0=osE('EPREFIX')
p1='/usr/share/sage/src'
p2=osE('MY_SAGE_IS_HERE')
include_0=find_sage_include_dir( p0, p1, p2 )
p4='/usr/include/csage/'
include_1=find_sage_include_dir( p0, p4, None )
ext_modules = \
 [
  Extension
   (
    "flint_sage", ["flint.pyx"], libraries=['flint','csage'], 
    include_dirs=[include_0,include_1]
   )
 ]

setup\
 (
  name = 'FLINT wrapper, accessed with Sage matrix or numbers',
  description = 'FLINT integer/rational matrice bindings',
  long_description = 'For now, some FLINT algorithms for integer linear '+
   'algebra exposed to Python',
  platforms = ['64-bit supported by sage, for details scan python.flint.sage.README for technical req'],
  cmdclass = {'build_ext':build_ext},
  version='20140309',
  ext_modules = ext_modules,
  author_email = 'quest: http://tiny.cc/DKryskov -> "Hello world"'+
   ' -> current address at the bottom',
  url='https://github.com/krisk0/razin',
  author='Крыськов Денис',
  license = 'GPL',
  classifiers=['Topic :: Scientific/Engineering :: Mathematics']
 )
