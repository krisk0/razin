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
IUSE="-static doc mpir"

DEPEND="dev-libs/mpfr
	dev-libs/ntl
 mpir? ( sci-libs/mpir )
	!mpir? ( dev-libs/gmp )
 doc? ( dev-texlive/texlive-latex )
	"
RDEPEND="${DEPEND}"

src_prepare() {
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

	# Fix QA notice about missing so name
 export EXTRA_SHARED_FLAGS="-Wl,-soname,${FLINT_LIB}"

 local with_mpir
 use mpir && with_mpir="--with-mpir=${EPREFIX}/usr" ||
             with_mpir="--with-gmp=${EPREFIX}/usr"

 echo ./configure \
		"$with_mpir" \
		--with-mpfr="${EPREFIX}"/usr \
		--with-ntl="${EPREFIX}"/usr \
		--prefix="${EPREFIX}"/usr \
		CC=$(tc-getCC) \
		CXX=$(tc-getCXX) \
		AR=$(tc-getAR)

	./configure \
		"$with_mpir" \
		--with-mpfr="${EPREFIX}"/usr \
		--with-ntl="${EPREFIX}"/usr \
		--prefix="${EPREFIX}"/usr \
		CC=$(tc-getCC) \
		CXX=$(tc-getCXX) \
		AR=$(tc-getAR) || die "configure failed"
}

src_compile(){
	emake verbose
	ln -s lib"${PN}$(get_libname ${PV})" lib"${PN}$(get_libname)"
	use doc && emake -C doc/latex manual
}

src_install(){
	default
	dosym ${FLINT_LIB} /usr/$(get_libdir)/lib${PN}$(get_libname)
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
