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
import os,sys,re
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

# perl/sed vile magic converted to Python
def re_sub( d, sou, tgt ):
 e=re.compile(sou)
 return e.sub( tgt, d )

def zap_subr( d, subr ):
 # the following expression eats too much :
 #  'static.*?'+subr+'.*?λ }λ'
 # so doing replacement 'static.*?'+subr -> subr in stupid way
 p0=d.find(subr)
 p1=d[:p0].rfind('static')
 d=d[:p1]+d[p0:]
 e=re.compile( subr+'.*?λ }λ' )
 return e.sub( '', d )

tail_warning='''
I would be greatly irritated if you tell me that this algorithm is wrong
 WITHOUT GIVING EXAMPLE OF INPUT DATA that make it fail

Report bugs via Github mechanism or e-mail

My e-mail is in my blog, detailed information on how to get it is close to tail
 of setup.py'''

def sed_and_perl__goodbye( oN, iN ):
 '''
 This subroutine filters C file, removing what is not needed in production
  code
 
 sed and perl on my computer fail to properly substitute regular expessions 
 so I re-implement the functionality of delete_debug_C_code.sh in Python
 '''
 unwanted=re.compile( 'MPLUS|MMUL' )
 with open( oN, 'wb' ) as o, open( iN, 'rb' ) as i:
  d,r=i.readlines(),''
  for s in d:
   if unwanted.search(s):
    continue
   r += s.replace('\n','λ')
  r = re_sub ( r, ' *#' ,'#' )
  r = re_sub ( r, '#if .*?#endifλ' ,'' )
  r = re_sub ( r, 'λλλstatic', 'λλstatic' )
  for s in 'nmod_mat_print','vec_print','DKryskov_nmod_easy_zl':
   r = zap_subr( r, s )
  for i in range(10):
   r=r.replace( 'λλλ', 'λλ' )
  o.write( r.replace('λ','\n') )
  o.write( '/*'+tail_warning+'\n*/' )

sed_and_perl__goodbye( 'nmod_mat_HNF.c', 'nmod_mat_HNF-debug.c' )

setup\
 (
  name = 'FLINT wrapper, accessed with Sage matrix or numbers',
  description = 'FLINT integer/rational matrice bindings',
  long_description = 'For now, some FLINT algorithms for integer linear '+
   'algebra exposed to Python',
  platforms = ['64-bit supported by Sage, for details scan python.flint.sa'+
               'ge.README for technical req'],
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
