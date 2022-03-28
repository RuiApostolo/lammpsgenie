########
# module containing commonly used data reading and analysis tools
# Created by Michael Doig
# ported to python3 by Rui Apóstolo
# import command - to be pasted in program code
import gzip
"""
from commondata import readAll, getNatoms, readTS, getTSrange, getAtomType
from commondata import getAtomData, getTS, getAt, getMol, getTotMass, getCOMts
from commondata import readcombmasses, readAllGzip
"""
###################################################################
"""
List of modules
1) readAll
2) getNatoms
3) readTS
4) getTSrange
5) getAtomType
6) getAtomData
7) getTS
8) getAt
9) getMol
10) getTotMass
11) getCOMTS
12) readcombmasses
"""
#############################################################


def readAll(filename):
    """
    Reads the entire file in line-by-line and returns list of lines.

    Parameters
    ----------
    filename : str
        name of file to be read.

    Returns
    -------
    lines : list
        A list containing the lines of the given file.
    """

    with open(filename, "r") as ifile:
        lines = ifile.readlines()
    return lines


def readAllGzip(filename):
    """
    Reads the entire gzipped file in line-by-line and returns list of
    lines.

    Parameters
    ----------
    filename (string): name of file to be read.

    Returns
    -------
    list: a list containing the lines of given file.
    """

    with gzip.open(filename, "r") as ifile:
        lines = ifile.readlines()
    return lines


def getNatoms(lines):
    """
    Gets the total number of atoms in dump file.

    Parameters
    ----------
    lines (list): read list of lines.

    Returns
    -------
    int:
    last modified: 13/03/14

    Notes
    -----
    Fails hard if lines don't contain the target string.
    """

    for i, x in enumerate(lines):
        if lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i+1])
            return natoms
