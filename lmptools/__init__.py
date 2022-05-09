# __init__.py
#  import sys
from lmptools.commondata_p3 import getDumpTSRange, readTS, getTrajTSRange
from lmptools.atoms import getNatoms, getAtomType, getAtomData, getAllAtomData, getAtomRange, getAtomsByType, getTotalMass, getCOM
from lmptools.readfiles import readAll, readAllGzip
from lmptools.mergedatafiles_p3 import readInputFile, readTopology, readTopologies, limitsTopology, limitsAllTopologies, absoluteLimitsTopologies, shiftTopology, shiftTopologies, mergeTopologies, writeTopology
from lmptools._version import __version__
from lmptools._name import _name

#  sys.path.append('.')

#  __ALL__ = []
# TODO understand the implications of these imports and whether they're needed
