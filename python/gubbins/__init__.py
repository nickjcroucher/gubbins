#!/usr/bin/env python
# encoding: utf-8

"""
Imports into the `gubbins` namespace all fundamental
classes and methods for instantiating objects in the
for usage by client code.
"""

import sys
import os
import pkg_resources

###############################################################################
## Populate the 'gubbins' namespace

from gubbins import common

###############################################################################
## PACKAGE METADATA

__project__ = "Gubbins"

try:
    try:
        __homedir__ = __path__[0]
    except AttributeError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except IndexError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
except OSError:
    __homedir__ = None
except:
    __homedir__ = None

__author__ = "Andrew Page, Nicholas Croucher, Aidan Delaney, Christoph Puethe and Simon Harris"
__copyright__ = "Copyright 2020 Wellcome Trust Sanger Institute and Imperial College London"
__license__ = """
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

def version():
    os.environ["PATH"] = os.environ["PATH"] + ":/usr/lib/gubbins/"
    program_version = ""
    try:
        program_version = str(pkg_resources.get_distribution(__project__).version)
    except pkg_resources.RequirementParseError:
        pass
    return "%s" % program_version

__version__ = version()

def description():
    return "%s %s" % (__project__, version())

if __name__ == "__main__":
    sys.stdout.write("%s\n" % description())

