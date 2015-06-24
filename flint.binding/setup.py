# -*- coding: utf-8

# This program is part of RAZIN

# Copyright Денис Крыськов 2014
# License: GNU General Public License (GPL)

# This program locates Sage library, crunches C code, then builds flint_sage 
#  Python module

# If this program fails to find sage include, set MY_SAGE_IS_HERE to point to 
#  directory containing sage/rings/integer

'''
Error message
 x86_64-pc-linux-gnu-gcc: error: : No such file or directory
is harmless. And I don't know where it comes from (it occurs before my code
in setup.py takes control)
'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os,sys,re,subprocess
osE=os.environ.get
write=sys.stdout.write
cout_please=subprocess.check_output
delete_us=[]

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

head_warning='This file is auto-generated from '+\
 'C/nmod_mat/nmod_mat_HNF-debug.cλ'
head_trigger='ast, it comes this wayλ'

def sed_and_perl__goodbye( oN, iN ):
 '''
 This subroutine filters C file, removing stuff not needed in production
  code
 
 I fought sed and perl. I lost the fight and say goodbye to them. Let us say
  RAZIN no longer depends on sed or perl caprice, and total code size decreased
 '''
 unwanted=re.compile( 'MPLUS|MMUL| ASSERT' )
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
  for s in 'nmod_mat_print','vec_print':
   r = zap_subr( r, s )
  for i in range(10):
   r=r.replace( 'λλλ', 'λλ' )
  if r[-1] != 'λ':
   r += 'λ'
  o.write( r.replace('λ','\n') )
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

def find_flint_so_in_dir(d):
 # user forgot to specify FLINT shared object name, however specified directory
 # Which means user wants any libflint.so in the directory
 import glob
 for x in glob.glob(d+'/libflint.so*'):
  return x
 print 'Warning: no file matching libflint.so* in '+d
 return 0

sed_and_perl__goodbye( 'nmod_mat_HNF.c', 'C/nmod_mat/nmod_mat_HNF-debug.c' )

'''
 sage-on-gentoo ebuilds put *.p* into $EPREFIX/usr/share/sage
'''

p0=osE('EPREFIX')
p1='/usr/share/sage/src'
p2=osE('MY_SAGE_IS_HERE')
include_0=find_sage_include_dir( p0, p1, p2 )
'''
With some versions of Sage .pyx fails to compile, error message: 
 can't find ccobject.h
Arranging extra include directory
'''
# /usr/share/sage/src -> /usr/share/sage/src/sage/ext/ccobject.h
find_ccobject_please=include_0+'/sage/ext/'
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
  flint_so=find_flint_so_in_dir( my_so )
  if 0 != flint_so:
   bad_path=0
 if bad_path:
  die( 'no such file or directory '+my_so+', or FLINT .so not found' )

print 'libraries taken via -l:',libraries
print 'libraries taken as extra objects:',extra_objects

def create__fmpz_mat_HNF(h_file,out_file):
 if os.system( "grep -q 'void fmpz_mat_hnf' '%s'" % h_file ):
  with open(out_file,'wb') as f:
   pass
 else:
  os.system( "cp '%s.in' '%s'" % (out_file,out_file) )
 global delete_us
 delete_us.append(out_file)

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

setup_automagic='# This file automagically generated by setup.py from %s\n\n'

def ctypedef( c ):
 try:
  d=cout_please( c, shell=True )
 except:
  print 'Failed to execute command\n%s\n' % c
  raise
 d=d.replace(';','').strip().replace('typedef',' ctypedef')
 return d

def mp_limb_t( gmpH, mpfrH ):
 ' invoke C pre-processor to extract 2 type definitions from gmp.h '
 c_pattern="%s -E %s|grep '%s;'|grep typedef"
 if not os.path.isfile( gmpH ):
  die('Failed to find gmp.h')
 global delete_us
 delete_us.append("mp_limb_t.pyx")
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

def slong( p ):
 if not os.path.isfile( p ):
  die('Failed to find FLINT includes')
 tU,tS=take_define( p, 'ulong' ),take_define( p, 'slong' )
 p=p.replace( '/flint.h', '/fmpz.h' )
 c="%s -E %s|grep 'fmpz;'|grep typedef" % (CC,p)
 fmpz=ctypedef( c )
 p = p[ :p.rfind('/') ] + '/*'
 global delete_us
 delete_us.append("slong.pyx")
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
 mp_limb_t( x+'/gmp.h', x+'/mpfr.h' )
 slong( x+'/flint/flint.h' )

def mullow_n( so ):
 global delete_us
 if 0==so:
  so='-lflint'
 if p0==None:
  incl=' '
  libr=' '
 else:
  incl=' "-I'+p0+'/usr/include"'
  libr=' "'+p0+'"'
 rc=os.system( CC+incl+' C/mullow_n.c -c -omullow_n.o -O2 -g' )
 if rc:
  die( 'mullow_n(): compilation problem. FLINT header not found?' )
 delete_us.append('mullow_n.o')
 exe='mullow.exe'
 rc=os.system( CC+' mullow_n.o '+so+libr+' -lgmp -o'+exe )
 if rc:
  # TODO: assume gmp version is wrong, extract gmp .so name from flint .so with
  #  objdump -p or some other way, re-link
  return 0
 delete_us.append(exe)
 rc=os.system('./%s %s' % (exe,(int(os.stat(exe).st_atime)|1) % 2**32))
 return 0==rc

def MulMod_2x( p ):
 if mullow_n( flint_so ):
  smart='smaart'
 else:
  smart='stupid'
 tgt=p+'_positive.c'
 os.system( "cp %s.%s %s" % (p,smart,tgt) )
 global delete_us
 delete_us.append(tgt)

CC=osE('CC','gcc')
include_dirs=[include_0,include_1,find_ccobject_please]
if p0 == None:
 slong__mp_limb_t( '/usr/include' )
else:
 slong__mp_limb_t( p0+'/usr/include' )
 include_dirs.append( p0+'/usr/include' )

MulMod_2x( 'C/fmpz/MulMod_2x' )

ext_modules = \
 [
  Extension
   (
    "flint_sage", ['flint.pyx'], 
    include_dirs=include_dirs,
    libraries=libraries, 
    extra_objects=extra_objects,
    library_dirs=library_dirs,
    runtime_library_dirs=runtime_library_dirs
   )
 ]

# tiny.cc not always works, good luck with email quest
setup\
 (
  name = 'RAZIN',
  description = 'Integer linear algebra library, built by Денис Крыськов '+
   'on top of FLINT',
  long_description = 'some FLINT and my home-brewed algorithms for '+
   'linear algebra and number theory exposed to Python',
  platforms = ['64-bit supported by Sage, look into python.flint.sage.README '+
                'for technical req'],
  cmdclass = {'build_ext':build_ext},
  version='RAZIN_version',
  ext_modules = ext_modules,
  author_email = 'quest: http://tiny.cc/DKryskov -> "Hello world"'+
   ' -> current address at the bottom',
  url='https://github.com/krisk0/razin',
  author='Крыськов Денис',
  license = 'GPL',
  classifiers=['Topic :: Scientific/Engineering :: Mathematics']
 )

if osE('ZAP_AUTOGENERATED'):
 write('deleting ')
 for x in delete_us:
  write(x+' ')
  os.remove(x)
 print
