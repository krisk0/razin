# -*- coding: utf-8

# This program is part of RAZIN

# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# This program locates Sage library, crunches C code, then builds flint_sage 
#  Python module

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
I would be greatly irritated if you tell me that the algorithms implemented in
 this file are wrong WITHOUT GIVING EXAMPLE OF INPUT DATA that make them fail

Report bugs via Github mechanism or e-mail

My e-mail is in my blog, detailed information on how to get it is close to tail
 of setup.py'''

head_warning='This file is auto-generated from C/nmod_mat_HNF-debug.c'+\
 'λλ#define NDEBUG 1'
head_trigger='ast, it comes this wayλ'

def sed_and_perl__goodbye( oN, iN ):
 '''
 This subroutine filters C file, removing what is not needed in production
  code
 
 I fought sed and perl. I lost the fight and say goodbye to them. Let us say
  RAZIN no longer depends on sed or perl caprice
 '''
 unwanted=re.compile( 'MPLUS|MMUL' )
 with open( oN, 'wb' ) as o, open( iN, 'rb' ) as i:
  d,r=i.readlines(),''
  for s in d:
   if unwanted.search(s):
    continue
   r += s.replace('\n','λ')
  r = re_sub ( r, ' *#' ,'#' )
  r = re_sub ( r, '#define BUG.*?λ', '' )
  r = re_sub ( r, '#if .*?#endifλ' , '' )
  r = r.replace( '#include "nmod', '#include "C/nmod_mat/nmod' )
  r = r.replace( head_trigger, head_trigger+'λ// '+head_warning+'λ' )
  #r = re_sub ( r, 'λλλstatic', 'λλstatic' )
  for s in 'nmod_mat_print','vec_print','DKryskov_nmod_easy_zl':
   r = zap_subr( r, s )
  for i in range(10):
   r=r.replace( 'λλλ', 'λλ' )
  if r[-1] != 'λ':
   r += 'λ'
  o.write( r.replace('λ','\n') )
  o.write( '\n#undef NDEBUG\n\n' )
  o.write( '/*'+tail_warning+'\n*/' )

sed_and_perl__goodbye( 'nmod_mat_HNF.c', 'C/nmod_mat/nmod_mat_HNF-debug.c' )

'''
 sage-on-gentoo ebuilds put *.p* into $EPREFIX/usr/share/sage
'''

p0=osE('EPREFIX')
p1='/usr/share/sage/src'
p2=osE('MY_SAGE_IS_HERE')
include_0=find_sage_include_dir( p0, p1, p2 ) # under Gentoo they are at 
p4='/usr/include/csage/' # Sage .h should be in $EPREFIX/usr/include/csage/
include_1=find_sage_include_dir( p0, p4, None )

# Link to flint dynamically or statically
# I failed to link statically under Linux, however leave code here 
# todo: arrange it so libflint.a and mpfr.a get statically linked
libraries=['flint','csage']
extra_objects,library_dirs,runtime_library_dirs=[],[],[]
my_so=osE('MY_FLINT_IS_HERE')
if my_so != None:
 bad_path=1
 if os.path.isfile(my_so):
  bad_path=0
  flint_so=my_so
  if flint_so[-2:] == '.a':
   extra_objects=[flint_so]
   libraries=[flint_so,'csage']# todo: maybe append mpfr.a to libraries
  else:
   libraries=[flint_so,'csage']
   # flint*.so will want to find mpfr*.so and gmp*.so shared libs,
   #  assume they are all in same directory. Add it to runtime_library_dirs
   runtime_library_dirs=[ os.path.dirname( flint_so ) ]
 if os.path.isdir(my_so):
  library_dirs=[ my_so ]
  runtime_library_dirs=[ my_so ]
  bad_path=0
 if bad_path:
  print 'no such file or directory',my_so
  sys.exit(1)

print 'shared libraries:',libraries
print 'static libraries:',extra_objects

ext_modules = \
 [
  Extension
   (
    "flint_sage", ["flint.pyx"], 
    include_dirs=[include_0,include_1],
    libraries=libraries, 
    extra_objects=extra_objects,
    library_dirs=library_dirs,
    runtime_library_dirs=runtime_library_dirs
   )
 ]

setup\
 (
  name = 'FLINT wrapper, accessed with Sage matrix or numbers',
  description = 'FLINT integer/rational matrice bindings',
  long_description = 'For now, some FLINT and my home-brewed algorithms for '+
   'integer linear algebra exposed to Python',
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
