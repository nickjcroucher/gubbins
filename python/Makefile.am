EXTRA_DIST=gubbins/*

all-local:
	pip install .

install-exec-local:
	pip install .

uninstall-local:
	rm -rf $(pythondir)/*gubbins*

clean-local:
	rm -rf $(pythondir)/build

check-local:
	pytest --cov-report xml --cov=gubbins gubbins/tests/
