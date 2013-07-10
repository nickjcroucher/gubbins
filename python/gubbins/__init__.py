#! /usr/bin/env python

"""
Imports into the `gubbins` namespace all fundamental
classes and methods for instantiating objects in the
for usage by client code.
"""

import sys
import os

###############################################################################
## Populate the 'gubbins' namespace

from gubbins import common


###############################################################################
## PACKAGE METADATA

__project__ = "Gubbins"
__version__ = "0.1"

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

__author__ = "Andrew J. Page, Nicholas Croucher, Aidan Delaney and Simon Harris"
__copyright__ = "Copyright 2013 Wellcome Trust Sanger Institutue"
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
PACKAGE_VERSION = __version__

def description():
    if __revision__.is_available:
        revision_text = " (%s)" % str(__revision__)
    else:
        revision_text = ""
    return "%s %s%s" % (__project__, __version__, revision_text)

if __name__ == "__main__":
    sys.stdout.write("%s\n" % description())

