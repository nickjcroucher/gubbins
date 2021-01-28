EXTRA_DIST=gubbins/*

all-local:
	${PYTHON} setup.py build


install-exec-local:
	${PYTHON} setup.py install

uninstall-local:
	rm -rf $(pythondir)/*gubbins*

clean-local:
	$(PYTHON) setup.py clean --all

check-local:
	$(PYTHON) setup.py test
