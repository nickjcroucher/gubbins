EXTRA_DIST=gubbins/* setup.py scripts/*

all-local:
	${PYTHON} setup.py build


install-exec-local:
	${PYTHON} setup.py install  --root=$(DESTDIR) --install-purelib=$(pythondir) --install-scripts=/usr/bin

uninstall-local:
	rm -rf $(pythondir)/*gubbins*

clean-local:
	$(PYTHON) setup.py clean --all

check-local:
	echo "Unimplemented"
