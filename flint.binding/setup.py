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
import os,sys,re,subprocess
osE=os.environ.get
write=sys.stdout.write
cout_please=subprocess.check_output

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

def die(x):
 print x
 sys.exit(1)

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
  RAZIN no longer depends on sed or perl caprice, and total code size decreased
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

def copy_pyx_file( oN, iN, tail ):
 with open( oN, 'wb' ) as o:
  o.write('include "%s"' % iN)
  o.write('\n\n' + tail + '\n')

def fix_include( c, ii ):
 ext_file_exist( c )
 p='''sed -e 's:#include "%s.h":#include "flint/%s.h":' -i %s'''
 for i in ii:
  os.system( p % (i,i,c) )

sed_and_perl__goodbye( 'nmod_mat_HNF.c', 'C/nmod_mat/nmod_mat_HNF-debug.c' )

'''
 sage-on-gentoo ebuilds put *.p* into $EPREFIX/usr/share/sage
'''

p0=osE('EPREFIX')
p1='/usr/share/sage/src'
p2=osE('MY_SAGE_IS_HERE')
include_0=find_sage_include_dir( p0, p1, p2 )
p4='/usr/include/csage/' # Sage .h should be in $EPREFIX/usr/include/csage/
include_1=find_sage_include_dir( p0, p4, None )

# Link to flint dynamically or statically
# I failed to link statically under Linux, however leave code here 
libraries=['flint','csage']
extra_objects,library_dirs,runtime_library_dirs=[],[],[]
flint_so=0
my_so=osE('MY_FLINT_IS_HERE')
if my_so != None:
 bad_path=1
 if os.path.isfile(my_so):
  bad_path=0
  flint_so=my_so
  extra_objects=[flint_so]
  runtime_library_dirs=[ os.path.dirname( flint_so ) ]
  libraries=['csage']
 if os.path.isdir(my_so):
  library_dirs=[ my_so ]
  runtime_library_dirs=[ my_so ]
  bad_path=0
 if bad_path:
  die( 'no such file or directory '+my_so )

print 'libraries taken via -l:',libraries
print 'libraries taken as extra objects:',extra_objects

def create__fmpz_mat_HNF(h_file,out_file):
 if os.system( "grep -q 'void fmpz_mat_hnf' '%s'" % h_file ):
  with open(out_file,'wb') as f:
   pass
 else:
  os.system( "cp '%s.in' '%s'" % (out_file,out_file) )

# create fmpz_mat_HNF.pyx: either empty or defining AlexBest_hnf()
fmpz_mat__h='/usr/include/flint/fmpz_mat.h'
if p0!=None:
 fmpz_mat__h=p0+fmpz_mat__h
create__fmpz_mat_HNF( fmpz_mat__h, 'fmpz_mat_HNF.pyx' )

def ext_file_exist(x):
 if not os.path.isfile(x):
  die( 'ext source file not found: '+x )

def define_extern(f, p, e):
 ext_file_exist(f)
 if os.system( "sed -e 's:%s:&\\n%s;\\n:' -i %s" % (p,e,f) ):
  die( 'failed to add definition %s to file %s' (e,f) )

def copy_func( r, s, f ):
 t=extract_func( s, f )
 with open(r,'ab') as f:
  f.write('\n'+t)

def extract_func( i, n ):
 ext_file_exist(i)
 m=''
 with open(i,'rb') as f:
  for x in f.readlines():
   m += x.rstrip()+'λ'    # zap trailling spaces
 p0=m.find(n)
 m=m[p0:]
 # head cut-off, now find tail
 p2=m.find( 'λ}λ' )
 assert p2>0
 return m[:p2+3].replace( 'λ', '\n' )

'''
I don't want RAZIN to contain unmodified code from another repository.

I don't want to replace official version of FLINT such as flint-2.4.4.tar.gz 
 with contents of some experimental repository (such as Alex Best's 
 https://github.com/alexjbest/flint2).

However I want to optionally wrap some experimental functions (such as Alex
 Best implementaion of Pernet-Stein double det HNF).
 
Dirty hack below optionally replaces flint.pyx with another file; removes the
 new file after compilation
''' 

FLINT_PYX="flint.pyx"
FLINT_PYX_copied=0
EXTRA_0="C/fmpz_mat/hnf_pernet_stein.c"
if os.path.isfile(EXTRA_0):
 EXTRA_1="C/fmpz_mat/hnf_modular.c"
 EXTRA_2="C/fmpz_mat/hnf_xgcd.c"
 EXTRA_3="C/fmpz_mat/hnf_classical.c"
 EXTRA_4="C/fmpz_mat/hnf_minors.c"
 # patch files EXTRA_?
 fix_include(EXTRA_0, ['fmpz_mat','fmpq_mat','perm'] )
 define_extern(EXTRA_0, 'perm.h"',
  'void @_hnf_modular(@_t,const @_t,const fmpz_t)'.replace('@','fmpz_mat') )
 define_extern(EXTRA_0, 'perm.h"',
  'void @_hnf_xgcd(@_t,const @_t)'.replace('@','fmpz_mat') )
 define_extern(EXTRA_0, 'perm.h"',
  'slong _nmod_mat_rref(nmod_mat_t A, slong* pivots_nonpivots)' )
 define_extern(EXTRA_0, 'perm.h"',
  'void @_hnf_classical(@_t H, const @_t A)'.replace('@','fmpz_mat') )
 define_extern(EXTRA_0, 'perm.h"',
  'void @_hnf_minors(@_t H, const @_t A)'.replace('@','fmpz_mat') )
 fix_include(EXTRA_1, ['fmpz_mat'] )
 fix_include(EXTRA_2, ['fmpz_mat'] )
 fix_include(EXTRA_3, ['fmpz_mat'] )
 fix_include(EXTRA_4, ['fmpz_mat'] )
 FLINT_PYX_copied=1
 FLINT_PYX="remove_me.pyx"
 copy_pyx_file( FLINT_PYX, "flint.pyx", 'include "wrap_HNF_pernet_stein.pyx"' )
 #ImportError: /path/to/flint_sage.so: undefined symbol: _nmod_mat_rref
 Ufunc='_nmod_mat_rref'
 if 0 == os.system( 'grep -q %s %s' % (Ufunc,EXTRA_0) ):
  copy_func( EXTRA_0, 'C/nmod_mat/rref.c', 'slongλ'+Ufunc )
  os.system( ("sed -e 's-@-U@-g' -i"+EXTRA_0).replace('@',Ufunc) )

setup_automagic='# This file automagically generated by setup.py from %s\n\n'

def ctypedef( c ):
 try:
  d=cout_please( c, shell=True )
 except:
  print 'Failed to execute command\n%s\n' % c
  raise
 d=d.replace(';','').strip().replace('typedef',' ctypedef')
 return d

def mp_limb_t( gmpH, mpfrH, CC ):
 ' invoke C pre-processor to extract 2 type definitions from gmp.h '
 c_pattern="%s -E %s|grep '%s;'|grep typedef"
 if not os.path.isfile( gmpH ):
  die('Failed to find gmp.h')
 with open("mp_limb_t.pyx",'wb') as o:
  o.write( setup_automagic % (gmpH+' and '+mpfrH) )
  o.write( "cdef extern from 'gmp.h':\n" )
  for t in 'mp_limb_t','mp_limb_signed_t','mp_bitcnt_t':
   o.write( ctypedef( c_pattern % (CC,gmpH,t) )+'\n' )
  o.write( "\ncdef extern from 'mpfr.h':\n" )
  for t in 'mpfr_prec_t', 'mpfr_sign_t', 'mpfr_exp_t', 'mpfr_uexp_t':
   o.write( ctypedef( c_pattern % (CC,mpfrH,t) )+'\n' )

def take_define( f, q ):
 with open(f,'rb') as i:
  for j in i:
   k=j.strip()
   while k.find('  ') >= 0:
    k=k.replace('  ',' ')
   k=k.split(' ')
   if len(k)==3 and k[0]=='#define' and k[1]==q:
    return k[2]
 die( 'failed to find definition of %s in %s' % (q,f) )

def slong( p, CC ):
 if not os.path.isfile( p ):
  die('Failed to find FLINT includes')
 tU,tS=take_define( p, 'ulong' ),take_define( p, 'slong' )
 p=p.replace( '/flint.h', '/fmpz.h' )
 c="%s -E %s|grep 'fmpz;'|grep typedef" % (CC,p)
 fmpz=ctypedef( c )
 p = p[ :p.rfind('/') ] + '/*'
 with open("slong.pyx",'wb') as o:
  o.write( setup_automagic % p)
  o.write( "cdef extern from 'flint/flint.h':\n" )
  o.write( ' ctypedef %s ulong\n' % tU )
  o.write( ' ctypedef %s slong\n' % tS )
  o.write( "\ncdef extern from 'flint/fmpz.h':\n" )
  o.write( fmpz+'\n' )

def slong__mp_limb_t( x ):
 '''
  assume FLINT, GMP and MPFR headers share same directory. 
  Extract basic number type definitions 
 '''
 CC=os.environ.get('CC','gcc')
 mp_limb_t( x+'/gmp.h', x+'/mpfr.h', CC )
 slong( x+'/flint/flint.h', CC )

include_dirs=[include_0,include_1]
if p0 == None:
 slong__mp_limb_t( '/usr/include' )
else:
 slong__mp_limb_t( p0+'/usr/include' )
 include_dirs.append( p0+'/usr/include' )

ext_modules = \
 [
  Extension
   (
    "flint_sage", [FLINT_PYX], 
    include_dirs=include_dirs,
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

if FLINT_PYX_copied:
 os.remove(FLINT_PYX)
