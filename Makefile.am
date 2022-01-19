EXTRA_DIST=debian/* tests/*.h tests/data/*
SUBDIRS=src release python
ACLOCAL_AMFLAGS=-I m4
CFLAGS := -std=gnu99
# include "check-code-coverage" rule
include $(top_srcdir)/aminclude_static.am

source: dist
	cd release && make source

release: dist
	cd release && make ubuntu-binary
	cd release && make source

osx-package:
	rm -rf ${PACKAGE}_${VERSION}_osx_10.8.tar.gz
	rm -rf ${HOME}/files.lst
	rm -rf ${HOME}/${PACKAGE}-${VERSION}-inst
	./configure --prefix /usr
	make
	make DESTDIR=${HOME}/${PACKAGE}-${VERSION}-inst install
	cd ${HOME}/${PACKAGE}-${VERSION}-inst && mv ${HOME}/${PACKAGE}-${VERSION}-inst/${HOME}/* .
	cd ${HOME}/${PACKAGE}-${VERSION}-inst && rm -rf ${HOME}/${PACKAGE}-${VERSION}-inst/Users
	cd ${HOME}/${PACKAGE}-${VERSION}-inst && find . -type f -print  > ${HOME}/files.lst
	cd ${HOME}/${PACKAGE}-${VERSION}-inst && tar zcvf ${PACKAGE}_${VERSION}_osx_10.8.tar.gz -T ${HOME}/files.lst
