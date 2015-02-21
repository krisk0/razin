# Copyright 1999-2014 Gentoo Foundation
# Copyright      2014 Крыськов Денис 
# Distributed under the terms of the GNU General Public License v2

# Based on flint-20140925.ebuild, changed to take version
#  d7b9299a2399c7fe4f0e50e59ca6926ccac79f9a dated 01 oct 2014 instead of a 25 Sept
#  version

EAPI=5

inherit eutils multilib toolchain-funcs

DESCRIPTION="Fast Library for Number Theory"
HOMEPAGE="http://www.flintlib.org/"
SRC_URI="http://www.flintlib.org/$P.tar.gz"

RESTRICT=mirror
LICENSE="GPL-2"
SLOT=0
KEYWORDS="amd64 x86 x86-macos x64-macos x86-linux amd64-linux"
IUSE="doc gc ntl static-libs test"

RDEPEND="dev-libs/gmp
	dev-libs/mpfr
	gc? ( dev-libs/boehm-gc )
	ntl? ( dev-libs/ntl )"
DEPEND="${RDEPEND}
	doc? (
		app-text/texlive-core
		dev-texlive/texlive-latex
		dev-texlive/texlive-latexextra
	)"

src_prepare() {
 # $PN-20141004.patch fixes 20 Aug error in flint.h
 for i in $PN-20141004 nmod_mat_print_pretty ; do
  epatch "$FILESDIR"/$i.patch || die "patch $i does not apply"
 done
}

src_configure() {
 # fix >>>QA Notice: The following shared libraries lack a SONAME<<<

 local FLINT_LIB=lib"$PN$(get_libname ${PV})"
 #local FLINT_LIB=lib"$PN".so.2.4.5
 local ESF="-Wl,-soname,$FLINT_LIB"
 local p="$EPREFIX/usr"
 
	# configure supports CFLAGS but not CXXFLAGS; Makefile supports both
 # No need to feed C*FLAGS to configure

 ./configure \
		--prefix="$p" \
		--with-gmp="$p" \
		--with-mpfr="$p" \
		$(usex ntl "--with-ntl=$p" "") \
		$(use_enable static-libs static) \
		$(usex gc "--with-gc=$p" "") \
		CC=$(tc-getCC) \
		CXX=$(tc-getCXX)\
		AR=$(tc-getAR) \
		|| die 'configure failed'

 # append soname voodoo to EXTRA_SHARED_FLAGS value, change shared lib name
 sed -e "s:^EXTRA_SHARED_FLAGS=.*:& $ESF:" \
     -e "s:^FLINT_LIB=.*:FLINT_LIB=$FLINT_LIB:" \
     -i Makefile
}

src_compile() {
	emake verbose
	use doc &&	emake -C doc/latex
}

src_test() {
 # With static-libs disabled, test suite checks installed version in /usr/lib
 #  which is pointless
 use static-libs || die "To perform test, you need to enable static-libs"
	emake AT= QUIET_CC= QUIET_CXX= QUIET_AR= check
}

src_install() {
	emake DESTDIR="${D}" LIBDIR="$(get_libdir)" install
	einstalldocs
	use doc && dodoc doc/latex/flint-manual.pdf
}
