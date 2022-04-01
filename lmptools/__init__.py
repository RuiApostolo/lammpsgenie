# __init__.py
import sys
from .commondata_p3 import getNatoms, getTSrange, getAtomType, readTS, getTS, getAtomData
from .readfiles import readAll, readAllGzip

sys.path.append('.')
