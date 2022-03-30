########
# module containing commonly used data reading and analysis tools
# Created by Michael Doig
# ported to python3 by Rui Apóstolo
# import command - to be pasted in program code
import gzip
import re
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
        Name of file to be read.

    Returns
    -------
    lines : list
        A list containing the lines of the given file.
    """

    with open(filename, "r") as ifile:
        lines = ifile.read().splitlines()
    return lines


def readAllGzip(filename):
    """
    Reads the entire gzipped file in line-by-line and returns list of
    lines.

    Parameters
    ----------
    filename : str
        Name of file to be read.

    Returns
    -------
    lines: list of str
        A list containing the lines of given file.
    """

    with gzip.open(filename, "r") as ifile:
        lines = ifile.readlines()
    return lines


def getNatoms(lines):
    """
    Gets the total number of atoms in dump file.

    Parameters
    ----------
    lines : list of str
        Reads list of lines from dumpfile.

    Returns
    -------
    natoms : int
        Number of atoms.

    Notes
    -----
    Fails hard if lines don't contain the target string.
    """

    for i, x in enumerate(lines):
        if lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i+1])
            return natoms


def getTSrange(lines):
    """
    Reads dump file to get list of timestep numbers.

    Parameters
    ----------
    lines : list of str
        Read list of lines.

    Returns
    -------
    tsrange : list of int
        List of timestep numbers.

    Notes
    -----
    Fails hard if lines don't contain the target string.
    """

    tsrange = []
    for i, _ in enumerate(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            tsrange.append(int(lines[i+1]))
    return tsrange


# TODO here downwards
def getAtomType(filename):
    """
    Reads data file and returns a dictionary of the atom types, retrieved
    from comment after `#` in each line in the `Pair Coeff` section of
    the datafile.

    Parameters
    ----------
    filename : str
        Name of the datafile to be read.

    Returns
    -------
    atomnames : dict of str
        A dictionary with the LAMMPS atom type numbers as keys, and the
        custom atom type names as values.
        {1: 'CPS', 2: 'OCB'}
    """

    p = re.compile(r"(\s*\d*\s*)(?P<name>[a-zA-Z]*\w*)(\s+\w*)")
    with open(filename, "r") as datafile:
        lines = datafile.read().splitlines()
    for lineindex, line in enumerate(lines):
        if "atom types" in line:
            natomtypes = int(line.split()[0])
            print(f"Found {natomtypes} atom types.")
            atomnames = {}
        if "Pair Coeffs" in line:
            for i in range(1, natomtypes + 1):
                matches = re.match(p, lines[lineindex + i + 1].split("#")[1])
                if matches.group('name') != '':
                    atomnames[i] = matches.group('name')
                else:
                    atomnames[i] = str(i)
            return atomnames
    return {}


def readTS(lines, natoms, tsnum, atomnames, header=9):
    """
    Reads in one timestep `tsnum`.

    TODO
    """

    traj = {}
    nlines = natoms + header
    firstline = ((tsnum-1)*nlines)
    lastline = firstline + nlines-1
    # print('firstline/lastline '+str(firstline)+' / '+str(lastline))
    # for i, x in enumerate(lines[firstline:lastline]):
    for i in range(firstline, lastline):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1])
            traj[ts] = {}
            traj[ts]["atom"] = {}
            traj[ts]["boxsize"] = {}
            traj[ts]["boxx"] = {}
            traj[ts]["boxy"] = {}
            traj[ts]["boxz"] = {}
        if "BOX BOUNDS" in lines[i]:
            boxx = list(map(float, lines[i+1].split()))
            boxy = list(map(float, lines[i+2].split()))
            boxz = list(map(float, lines[i+3].split()))
            header = 9
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            traj[ts]["boxx"] = (boxx)
            traj[ts]["boxy"] = (boxy)
            traj[ts]["boxz"] = (boxz)
            # traj[ts]["boxsize"]['lx'] = lx
            # traj[ts]["boxsize"]['ly'] = ly
            # traj[ts]["boxsize"]['lz'] = lz
            traj[ts]["boxsize"] = (lx, ly, lz)
        if "xy xz yz pp pp" in lines[i]:
            boxx = list(map(float, lines[i+1].split()))
            boxy = list(map(float, lines[i+2].split()))
            boxz = list(map(float, lines[i+3].split()))
            header = 9
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            traj[ts]["boxx"] = (boxx)
            traj[ts]["boxy"] = (boxy)
            traj[ts]["boxz"] = (boxz)
            traj[ts]["boxsize"] = (lx, ly, lz)
        if lines[i].startswith("ITEM: ATOMS"):
            labels = lines[i].split()
            for j in range(natoms):
                atomval = lines[i+j+1].split()
                atomid = int(atomval[0])
                traj[ts]["atom"][atomid] = {}
                for label_idx, k in enumerate(labels[3:], start=1):
                    if k == 'type':
                        traj[ts]["atom"][atomid][k] = int(atomval[label_idx])
                    else:
                        traj[ts]["atom"][atomid][k] = float(atomval[label_idx])
                # Change from numerical atom type (from dump) to atom type
                # name (from data file)
                traj[ts]["atom"][atomid]["type"] = \
                    atomnames[traj[ts]["atom"][atomid]["type"]]
    return traj


def getTS(traj, ts1, ts2):
    """
    Select a partial timestep range from a trajectory.

    Parameters
    ----------
    traj : dict
        A trajectory 'object' from readTS. Keys are the timestep labels.
    ts1 : int or 'first'
        The timestep label for the start of the range.
        If 'first' uses the first timestep of the trajectory, regardless
        of number.
    ts2 : int or 'last'
        The timestep label for the end of the range.
        If 'last' uses the last timestep of the trajectory, regardless
        of number.

    Returns
    -------
    tsrange : list
        List of timestep labels.
    """

    # TODO: remove sorted and test if still functional
    if ts1 == 'first':
        ts1 = sorted(traj.keys())[0]
    if ts2 == 'last':
        ts2 = sorted(traj.keys())[-1]
    sortedts = sorted(traj.keys())
    tsrange = []
    for timestep in sortedts:
        if int(ts1) <= int(timestep) <= int(ts2):
            tsrange.append(timestep)
    return tsrange
