#/bin/sh
###############################################################

all:
	cd src; $(MAKE) -j 4

allsrc:
	cd src; $(MAKE)

clean:
	cd src; $(RM) .depends ; $(MAKE) clean
	cd tests; $(MAKE) clean

depends:
	cd src; $(MAKE) depends

install:
	cd src; $(MAKE) install

rebaseline:
	cd tests; $(MAKE) rebaseline

test:
	cd tests; $(MAKE) test


# end of file

