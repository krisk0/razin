# Copyright 1999-2013 Gentoo Foundation
# Copyright 2014 Денис Крыськов
# Distributed under the terms of the GNU General Public License v2
#  changes since mpir-2.6.0-r2:
#   get rid of patches since none of them apply
#   get rid of yasm, by removing the directory and patching Makefile.am
#   workaround ABI shell var problem

# This ebuild appears to be working, and test passes

# Known issues: for test to work, package need to be installed first:
#  emerge =mpir-2.7.0 && FEATURES=test emerge =mpir-2.7.0

EAPI=5

inherit autotools-utils eutils toolchain-funcs

DESCRIPTION="Library for arbitrary precision integer arithmetic (fork of gmp)"
HOMEPAGE="http://www.mpir.org/"
SRC_URI="http://www.mpir.org/${P}-alpha4.tar.bz2"

LICENSE="LGPL-3"
SLOT="0"
KEYWORDS="amd64 amd64-linux"
IUSE="+cxx cpudetection static-libs"

DEPEND="dev-lang/yasm"
RDEPEND=""

src_prepare() {
 [ -d yasm ] || die "yasm not found"
 rm -rf yasm || die "rm command failed, goodbye"
 sed -e 's:SUBDIRS += yasm::' -i Makefile.am || die 'patching Makefile failed'
	tc-export CC

	# In the same way there was QA regarding executable stacks
	# with GMP we have some here as well. We cannot apply the
	# GMP solution as yasm is used, at least on x86/amd64.
	# Furthermore we are able to patch config.ac.
	ebegin "Patching assembler files to remove executable sections"
	local i
	for i in $(find . -type f -name '*.asm') ; do
		cat >> $i <<-EOF
			#if defined(__linux__) && defined(__ELF__)
			.section .note.GNU-stack,"",%progbits
			#endif
		EOF
	done

	for i in $(find . -type f -name '*.as') ; do
		cat >> $i <<-EOF
			%ifidn __OUTPUT_FORMAT__,elf
			section .note.GNU-stack noalloc noexec nowrite progbits
			%endif
		EOF
	done
	eend
	eautoreconf
}

src_configure() {
	myeconfargs+=(
		$(use_enable cxx)
		$(use_enable cpudetection fat)
	)

 # ABI problem: Gentoo sets ABI=amd64 which confuses GMP and MPIR.
 # Workaround taken from gmp-6.0.0a ebuild
 export GMPABI=64
 mv configure configure.wrapped || die
	cat <<-\EOF > configure
	#!/bin/sh
	exec env ABI=$GMPABI "$0.wrapped" --with-system-yasm  --enable-gmpcompat "$@"
	EOF
	chmod a+rx configure
 
	autotools-utils_src_configure
}
