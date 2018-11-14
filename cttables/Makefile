#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2008/09/15 23:31:40 $
#  $Revision: 1.6 $
#
###############################################################
build_python=0

all:
	cd src; gmake

allsrc:
	cd src; gmake

clean:
	cd src; $(RM) .depends ; gmake clean
	cd tests; gmake clean

depends:
	cd src; gmake depends

install:
	cd src; gmake install

rebaseline:
	cd test; gmake rebaseline

test:
	cd tests; gmake test

# end of file

