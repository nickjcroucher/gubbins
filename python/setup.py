#!/usr/bin/env python3

import setuptools

setuptools.setup(
    name='gubbins',
    version = open('VERSION').read().strip(),
    description='Frontend to the Gubbins BioInformatics tool',
    author='Andrew Page, Nicholas Croucher, Aidan Delaney, Christoph Puethe and Simon Harris',
    author_email='n.croucher@imperial.ac.uk',
    url='https://github.com/sanger-pathogens/gubbins/',
    packages=setuptools.find_packages(),
    entry_points={
        "console_scripts": [
            "run_gubbins.py = gubbins.run_gubbins:main",
        ]
    },
    scripts=[
        'scripts/generate_ska_alignment.py',
        'scripts/extract_gubbins_clade.py',
        'scripts/mask_gubbins_aln.py',
        'scripts/gubbins_alignment_checker.py'
    ],
    tests_require=[
        "pytest >= 4.6",
        "wheel >= 0.34",
        "biopython >= 1.59",
        "dendropy  >= 4.0.2",
        "multiprocess >= 0.70",
        "scipy >= 1.5.3",
        "numpy >= 1.19",
        "ska >= 1.0"
    ],
    long_description="""\
      Gubbins is a tool that generates a reconstruction of
      a microbial strain's recent evolutionary history through
      identifying imports of divergent sequence through recombination
      and generating a phylogeny from the remaining clonal frame.
      """,
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
        "biopython >= 1.59",
        "dendropy  >= 4.0.2",
        "multiprocess >= 0.70",
        "scipy >= 1.5.3",
        "numpy >= 1.19"
    ],
    license="GPL"
)
