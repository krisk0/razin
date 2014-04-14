# Copyright 1999-2014 Gentoo Foundation
# Copyright      2014 Крыськов Денис 
# Distributed under the terms of the GNU General Public License v2

# This file is a fork of flint-2.4.3-r1.ebuild found at RAZIN git 
#  (https://github.com/krisk0/razin)

# Main difference from flint-2.4.3.ebuild: two extensions (ANTIC and BLAND) are
#  built when ext flag is set

# Ups, this ebuild does not build any extension due to a bug in FLINT 
#  configure/make voodoo

EAPI=5

inherit eutils multilib flag-o-matic toolchain-funcs

DESCRIPTION="Fast Library for Number Theory"
HOMEPAGE="http://www.flintlib.org/"

RESTRICT=mirror
LICENSE="GPL-2"
SLOT=0
KEYWORDS="amd64 x86 amd64-linux x86-linux"
IUSE="-static shared doc mpir ntl static-pic ext"
auth0=fredrik-johansson
auth1=wbhart
SRC_URI="
http://www.flintlib.org/$P.tar.gz
https://github.com/$auth0/bland/archive/master.zip -> $auth0-bland.zip
https://github.com/$auth1/antic/archive/trunk.zip -> $auth1-antic.zip
"

DEPEND="dev-libs/mpfr
 ntl? ( dev-libs/ntl )
 mpir? ( sci-libs/mpir )
 !mpir? ( dev-libs/gmp )
 doc? ( dev-texlive/texlive-latex )
 "
RDEPEND="$DEPEND"

src_prepare() {
 use static || use shared ||  
  {
   einfo 'It makes no sense to install flint without static or shared library'
   einfo 'Set static USE flag if you want static library'
   einfo 'Set shared USE flag if you want shared library'
   die 'Change USE flags'
  }
 sed -e "s:\${GMP_DIR}/lib\":\${GMP_DIR}/$(get_libdir)\":" \
  -e "s:\${MPFR_DIR}/lib\":\${MPFR_DIR}/$(get_libdir)\":" \
  -e "s:\${NTL_DIR}/lib\":\${NTL_DIR}/$(get_libdir)\":" \
  -i configure
 sed -i "s:\$(DESTDIR)\$(PREFIX)/lib:\$(DESTDIR)\$(PREFIX)/$(get_libdir):g" \
  Makefile.in
}

src_configure() {
 # handwritten script, needs extra stabbing
 export FLINT_LIB=lib"${PN}$(get_libname ${PV})"
 if ! (use static) ; then
  sed -i "s:STATIC=1:STATIC=0:" configure
 fi

 # Fix QA notice about missing so name; add rpath so flint .so will find its 
 #  dependencies later
 export EXTRA_SHARED_FLAGS="-Wl,-soname,${FLINT_LIB} '-Wl,-rpath,"
 EXTRA_SHARED_FLAGS="$EXTRA_SHARED_FLAGS""$EPREFIX"/usr/lib"'"
 # this should work even if $EPREFIX contains spaces. However flint build system
 #  is consisitently broken with respect to space inside filename
 
 # add -fPIC to CFLAGS if building position-indepependent .a
 use static-pic && CFLAGS="$CFLAGS -fPIC"
 
 local with_mpir with_ntl with_mp with_ext
 use mpir && with_mpir="--with-mpir=$EPREFIX/usr" ||
             with_mpir="--with-gmp=$EPREFIX/usr"

 use ntl && with_ntl="--with-ntl=$EPREFIX/usr"
 
 use ext && 
  with_ext="--extensions=\"$WORKDIR/bland-master $WORKDIR/antic-trunc\""

 # configure dislikes "$empty_string", two potentially empty strings are here:
 #  $with_ntl and $with_ext.
 # $with_ntl is moved to bottom, and branch saves from empty $with_ext
 use ext && 
  {
   ./configure \
    "$with_mpir"\
    "$with_ext"\
    --with-mpfr="$EPREFIX/usr"\
    --prefix="$EPREFIX/usr"\
    CC=$(tc-getCC)\
    CXX=$(tc-getCXX)\
    AR=$(tc-getAR)\
    "$with_ntl"\
    || die "configure failed"
  }
 use ext ||
  {
   ./configure \
    "$with_mpir"\
    --with-mpfr="$EPREFIX/usr"\
    --prefix="$EPREFIX/usr"\
    CC=$(tc-getCC)\
    CXX=$(tc-getCXX)\
    AR=$(tc-getAR)\
    "$with_ntl"\
    || die "configure failed"
  }
}

find_lib(){
 r=`ls "${EPREFIX}/usr/lib"|egrep lib$1.so\..*\..*\..*|head -1`
 [ -n "$r" ] && { echo "${EPREFIX}/usr/lib/$r"; return; }
 r=`ls "${EPREFIX}/usr/lib"|fgrep lib$1.so.|head -1`
 [ -n "$r" ] && { echo "${EPREFIX}/usr/lib/$r"; return; }
 r=`ls "${EPREFIX}/usr/lib"|egrep lib$1.so|head -1`
 [ -n "$r" ] && { echo "${EPREFIX}/usr/lib/$r"; return; }
 r=`ls "${EPREFIX}/usr/lib"|egrep lib$1.a*|head -1`
 [ -n "$r" ] && { echo "${EPREFIX}/usr/lib/$r"; return; }
 die "failed to find library $1 in ${EPREFIX}/usr/lib"
}

src_compile(){

 sed -e s:AT=@:AT=: -i Makefile

 # stupid linker takes mpfr library from /usr/lib64 if it is present there, 
 #  disrespecting -L${EPREFIX}/usr and LD_LIBRARY_PATH; -nodefaultlibs does 
 #  not work either
 # We want our mpfr from ${EPREFIX}/usr.

 # We do the linker's work and inject full path of the library we want
 sed -e s:-lmpfr:`find_lib mpfr`:g -i Makefile
 sed -e s:-lgmp:`find_lib gmp`:g -i Makefile
 use ntl && sed -e s:-lntl:`find_lib ntl`:g -i Makefile
 
 emake verbose
 ln -s lib"${PN}$(get_libname ${PV})" lib"${PN}$(get_libname)"
 use doc && emake -C doc/latex manual
}

src_install(){
 ! use shared && rm libflint*.so*
 use static && dolib.a libflint*.a
 default
 use shared && dosym ${FLINT_LIB} /usr/$(get_libdir)/lib${PN}$(get_libname)
 use doc &&
  {
   insinto /usr/share/doc/$PF
   doins doc/latex/*pdf || die 'failed to install pdf'
  }
}

src_test(){
 # build/interfaces/NTL-interface.o required for test
 emake build/interfaces/NTL-interface.o
 emake check AT="" QUIET_CC="" QUIET_CXX="" QUIET_AR=""
}
