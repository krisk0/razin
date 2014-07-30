EAPI=5

PYTHON_COMPAT=( python2_7 ) # ?

inherit distutils-r1

DESCRIPTION="Python wrapper for part of eigen"
HOMEPAGE="https://pypi.python.org/pypi/$PN"
MD5=883c251c6a85fd78a32a76c13d2e9828
fn=$P-1.tar.gz
SRC_URI="https://pypi.python.org/packages/source/m/$PN/$fn#md5=$MD5 -> $fn"

LICENSE=LGPL-3
SLOT=0

KEYWORDS="amd64"
RESTRICT=mirror
#TODO: specify DEPEND and RDEPEND
#DEPEND="dev-cpp/eigen ?boost

PATCHES=$FILESDIR/$PN-boost.path.patch

S=$WORKDIR/$P-1
