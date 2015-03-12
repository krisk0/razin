# Copyright 2014 Денис Крыськов
# License: GNU General Public License (GPL)

# Error message
#  x86_64-pc-linux-gnu-gcc: error: : No such file or directory
# is harmless and should be ignored

EAPI=5

PYTHON_COMPAT=( python2_7 ) # I only tested 2.7

inherit distutils-r1

DESCRIPTION="Python wrapper for FLINT, sage library required"
MYn=razin
HOMEPAGE=https://github.com/krisk0/$MYn
MYnv=$MYn-$PV-$PR
SRC_URI="$HOMEPAGE/archive/$PV-$PR.zip -> $MYnv.zip"

LICENSE=GPL
SLOT=0
KEYWORDS=amd64
RESTRICT=mirror

DEPEND=sci-mathematics/sage
RDEPEND="
         $DEPEND
         sci-mathematics/sage-clib
         sci-mathematics/flint
        "

S="$WORKDIR/$MYnv/flint.binding"
DOCS=( python.sage.flint.README )

src_prepare()
 {
  sed -e s:RAZIN_version:$PV-$PR: -i setup.py
 }
