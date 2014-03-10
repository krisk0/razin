# Copyright 2014 Крыськов Денис
# License: GNU General Public License (GPL)

EAPI=4

PYTHON_COMPAT=( python{2_7,3_1,3_2} ) # I only tested 2.7

inherit distutils-r1

DESCRIPTION="Python wrapper for FLINT, sage library required"
MYn=razin
HOMEPAGE=https://github.com/krisk0/$MYn
MYnv=${MYn}-$PV
SRC_URI="$HOMEPAGE/archive/$PV.zip -> $MYnv.zip"

LICENSE=GPL
SLOT=0
KEYWORDS="amd64"
RESTRICT=mirror

DEPEND=sci-mathematics/sage
RDEPEND="
         $DEPEND
         sci-mathematics/sage-clib
         sci-mathematics/flint
        "

S="$WORKDIR/$MYnv/flint.binding"
DOCS=( python.flint.sage.README )
