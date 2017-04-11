# Copyright 2014-2017 Денис Крыськов
# License: GNU General Public License (GPL)

EAPI=5

PYTHON_COMPAT=( python2_7 ) # I only tested 2.7

inherit distutils-r1

DESCRIPTION="Python wrapper for RAZIN, sage library required"
MYn=razin
HOMEPAGE=https://github.com/krisk0/$MYn
MYnv=$MYn-$PV-$PR
SRC_URI="$HOMEPAGE/archive/$PV-$PR.tar.gz -> $MYnv.tar.gz"

LICENSE=GPL
SLOT=0
KEYWORDS=amd64
RESTRICT=mirror

DEPEND=sci-mathematics/sage
RDEPEND=$DEPEND

S="$WORKDIR/$MYnv/flint.binding"
DOCS=( python.flint.sage.README )

src_prepare()
 {
  sed -e s:RAZIN_version:$PV-$PR: -i setup.py
 }
