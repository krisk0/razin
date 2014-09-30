# Copyright 1999-2014 Gentoo Foundation
# Distributed under the terms of the GNU General Public License v2
# $Header: /var/cvsroot/gentoo-x86/sci-mathematics/flint/flint-2.4.4-r1.ebuild,v 1.3 2014/08/14 16:22:26 phajdan.jr Exp $

# Hacked by Kryskov Denis to take commit a980f9695a8b03e7e1ebdf94e2f8630aec6a59de
#  dated 25 sept 2014 instead of official 2.4.4

EAPI=5

inherit eutils multilib toolchain-funcs

DESCRIPTION="Fast Library for Number Theory"
HOMEPAGE="http://www.flintlib.org/"
#SRC_URI="http://www.flintlib.org/${P}.tar.gz"
SHA='a980f9695a8b03e7e1ebdf94e2f8630aec6a59de'
SRC_URI="https://github.com/wbhart/${PN}2/archive/$SHA.zip -> 
 $P.zip"
S="$WORKDIR/${PN}2-$SHA"

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
 # 	${PN}-2.4.3-whitespaces.patch and ntl62.patch do not apply, skipped
	epatch "${FILESDIR}"/${PN}-2.4.3-libdir.patch \
		"${FILESDIR}"/$PN-2.4.3-cflags-ldflags.patch \
		"${FILESDIR}"/$PN-2.4.4-test.patch \
		"${FILESDIR}"/$P-latex.patch 
}

src_configure() {
 # added -fpermissive to CXXFLAGS, for NTL-interface.lo to compile

	./configure \
		--prefix="${EPREFIX}/usr" \
		--with-gmp="${EPREFIX}/usr" \
		--with-mpfr="${EPREFIX}/usr" \
		$(usex ntl "--with-ntl=${EPREFIX}/usr" "") \
		$(use_enable static-libs static) \
		$(usex gc "--with-gc=${EPREFIX}/usr" "") \
		CC=$(tc-getCC) \
		CXX="$(tc-getCXX) -fpermissive"\
		AR=$(tc-getAR) \
		|| die
}

src_compile() {
	emake verbose

	use doc &&	emake -C doc/latex
}

src_test() {
	emake AT= QUIET_CC= QUIET_CXX= QUIET_AR= check
}

src_install() {
	emake DESTDIR="${D}" LIBDIR="$(get_libdir)" install
	einstalldocs
	use doc && dodoc doc/latex/flint-manual.pdf
}
