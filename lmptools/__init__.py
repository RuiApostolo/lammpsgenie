# __init__.py
import sys
from .commondata_p3 import getDumpTSRange, readTS, getTrajTSRange
from .atoms import getNatoms, getAtomType, getAtomData
from .readfiles import readAll, readAllGzip

sys.path.append('.')
