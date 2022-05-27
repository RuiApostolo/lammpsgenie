# __init__.py
"""
This package is a collection of tools for handling LAMMPS data and dump files.

Please refer to more information in the README.md provided at:
github.com/RuiApostolo/lmpdtmrg
"""

from lammpsgenie.commondata import *
from lammpsgenie.atoms import *
from lammpsgenie.readfiles import *
from lammpsgenie.mergedatafiles import *
from lammpsgenie._version import __version__
from lammpsgenie._name import _name
