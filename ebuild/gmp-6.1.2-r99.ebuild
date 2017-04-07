# Copyright 1999-2017 Gentoo Foundation
#                2017 Денис Крыськов
# Un-cripple GMP by allowing it to optimize for microarch
# Distributed under the terms of the GNU General Public License v2

EAPI=5

inherit flag-o-matic eutils libtool multilib-minimal

MY_PV=${PV/_p*}
MY_PV=${MY_PV/_/-}
MY_P=${PN}-${MY_PV}
PLEVEL=${PV/*p}
DESCRIPTION="Library for arithmetic on arbitrary precision integers, rational numbers, and floating-point numbers"
HOMEPAGE="http://gmplib.org/"
SRC_URI="ftp://ftp.gmplib.org/pub/${MY_P}/${MY_P}.tar.xz
 mirror://gnu/${PN}/${MY_P}.tar.xz
 doc? ( http://gmplib.org/${PN}-man-${MY_PV}.pdf )"

LICENSE="|| ( LGPL-3+ GPL-2+ )"
# The subslot reflects the C & C++ SONAMEs.
SLOT="0/10.4"
KEYWORDS=amd64
IUSE="+asm doc cxx pgo static-libs"

DEPEND="sys-devel/m4
 app-arch/xz-utils"
RDEPEND=""

S=${WORKDIR}/${MY_P%a}

DOCS=( AUTHORS ChangeLog NEWS README doc/configuration doc/isa_abi_headache )
HTML_DOCS=( doc )
MULTILIB_WRAPPED_HEADERS=( /usr/include/gmp.h )

src_prepare() {
 [[ -d ${FILESDIR}/${PV} ]] && EPATCH_SUFFIX="diff" EPATCH_FORCE="yes" epatch "${FILESDIR}"/${PV}

 # note: we cannot run autotools here as gcc depends on this package
 elibtoolize

 epatch "${FILESDIR}"/${PN}-6.1.0-noexecstack-detect.patch

 # GMP uses the "ABI" env var during configure as does Gentoo (econf).
 # So, to avoid patching the source constantly, wrap things up.
 mv configure configure.wrapped || die
 printf '#!/usr/bin/env sh\nexec env ABI="${GMPABI}" "$0.wrapped" "$@"' > \
  configure
 # Patches to original configure might have lost the +x bit.
 chmod a+rx configure{,.wrapped}

 # multilib_src_configure() clobbers config.guess, so we run it here
 build_alias=`/bin/sh $S/config.guess` || die "failed to run config.guess"
 [ -z $build_alias ] && die "empty result from config.guess"
 einfo "guessed processor type: $build_alias"
}

multilib_src_configure() {
 # Because of our 32-bit userland, 1.0 is the only HPPA ABI that works
 # http://gmplib.org/manual/ABI-and-ISA.html#ABI-and-ISA (bug #344613)
 if [[ ${CHOST} == hppa2.0-* ]] ; then
  GMPABI="1.0"
 fi

 # ABI mappings (needs all architectures supported)
 case ${ABI} in
  32|x86)       GMPABI=32;;
  64|amd64|n64) GMPABI=64;;
  [onx]32)      GMPABI=${ABI};;
 esac
 export GMPABI

 #367719
 if [[ ${CHOST} == *-mint* ]]; then
  filter-flags -O?
 fi

 tc-export CC
 export ac_cv_host=$build_alias
 [ -z $ac_cv_host ] && die 'problem with build_alias'
 export ac_build_alias=$ac_cv_host
 ECONF_SOURCE="${S}" econf \
  --localstatedir="${EPREFIX}"/var/state/gmp \
  --enable-shared \
  $(use_enable asm assembly) \
  $(use_enable cxx) \
  $(use_enable static-libs static)
 unset ac_cv_host ac_build_alias build_alias
}

multilib_src_compile() {
 emake

 if use pgo ; then
  emake -j1 -C tune tuneup
  ebegin "Trying to generate tuned data"
  ./tune/tuneup | tee gmp.mparam.h.new
  if eend $(( 0 + ${PIPESTATUS[*]/#/+} )) ; then
   mv gmp.mparam.h.new gmp-mparam.h || die
   emake clean
   emake
  fi
 fi
}

multilib_src_test() {
 emake check
}

multilib_src_install() {
 emake DESTDIR="${D}" install

 # should be a standalone lib
 rm -f "${ED}"/usr/$(get_libdir)/libgmp.la
 # this requires libgmp
 local la="${ED}/usr/$(get_libdir)/libgmpxx.la"
 use static-libs \
  && sed -i 's:/[^ ]*/libgmp.la:-lgmp:' "${la}" \
  || rm -f "${la}"
}

multilib_src_install_all() {
 einstalldocs
 use doc && cp "${DISTDIR}"/gmp-man-${MY_PV}.pdf "${ED}"/usr/share/doc/${PF}/
}
