# Copyright 1999-2014 Gentoo Foundation
# Copyright      2014 Крыськов Денис 
# Distributed under the terms of the GNU General Public License v2

EAPI=5

inherit eutils multilib flag-o-matic toolchain-funcs

DESCRIPTION="Fast Library for Number Theory"
HOMEPAGE="http://www.flintlib.org/"
SRC_URI="http://www.flintlib.org/${P}.tar.gz"

RESTRICT=mirror
LICENSE="GPL-2"
SLOT=0
KEYWORDS="amd64 x86 amd64-linux x86-linux"
IUSE="-static shared doc mpir ntl static-pic"  
# static and static-pic: make sure .a is position-independent

DEPEND="dev-libs/mpfr
 ntl? ( dev-libs/ntl )
 mpir? ( sci-libs/mpir )
 !mpir? ( dev-libs/gmp )
 doc? ( dev-texlive/texlive-latex )
 "
RDEPEND="${DEPEND}"

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
 for p in nmod_mat_set_mod Abhinav_Baid__gcd nmod_mat_print_pretty \
   fmpz_txt_typo ; do
  epatch "$FILESDIR/$p.patch" || die "$p patch failed"
 done
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
 # this should work even if $EPREFIX contains spaces
 
 # add -fPIC to CFLAGS if building position-indepependent .a
 use static-pic && CFLAGS="$CFLAGS -fPIC"
 
 # add -nodefaultlibs so mpfr and other libraries will be taken from EPREFIX 
 #  rather than from /usr/lib64
 #EXTRA_SHARED_FLAGS="$EXTRA_SHARED_FLAGS -nodefaultlibs"
 # this trick failed, commented out

 local with_mpir with_ntl with_mp
 use mpir && with_mpir="--with-mpir=${EPREFIX}/usr" ||
             with_mpir="--with-gmp=${EPREFIX}/usr"

 use ntl && with_ntl="--with-ntl=${EPREFIX}/usr" || ntl=''

 # with empty $with_ntl in-between other parameters configure malfunctions,
 #  so $with_ntl moved to the bottom
 ./configure \
  "$with_mpir" \
  --with-mpfr="${EPREFIX}"/usr \
  --prefix="${EPREFIX}"/usr \
  CC=$(tc-getCC) \
  CXX=$(tc-getCXX) \
  AR=$(tc-getAR) \
  "$with_ntl"    \
  || die "configure failed"
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
