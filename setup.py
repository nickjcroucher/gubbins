#!/usr/bin/env python

from setuptools import setup

setup(
    name='gubbins',
    version='0.1',
    description='Frontend to the Gubbins BioInformatics tool',
    author='Andrew J. Page',
    author_email='ap13@sanger.ac.uk',
    url='https://github.com/andrewjpage/gubbins/',
    scripts=['run_gubbins.py'],
    long_description="""\
      Gubbins is a tool for tool for BioInformaticians that takes in a multi
      fasta alignment and detects recombination regions.  This package provides
      a simple front end to the Gubbins tool.
      """,
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    install_requires=[
        'Biopython >= 1.59',
        'Dendropy  >= 3.11.1',
        'Reportlab >= 2.5',
        ],
    license='GPL',
    )
