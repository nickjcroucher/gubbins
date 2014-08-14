EXTRA_DIST=debian/* tests/*.h tests/data/*
SUBDIRS=src release python
ACLOCAL_AMFLAGS=-I m4

source: dist
	cd release && make source
