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
    atomnames : dict
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


def readTS(lines, tsnum, atomnames, header=9):
    """
    Reads dumpfile and returns the atomic information of a timestep.

    Parameters
    ----------
    lines : list of str
        Reads list of lines from dumpfile.

    tsnum : int
        Number of the frame to read.

    atomnames : dict
        Dictionary of atom types and names. Use getAtomType().

    header : int, default=9
        Number of lines in the dumpfile 'header', from 'ITEM: TIMESTEP' to
        'ITEM: ATOM', inclusive.
    Returns
    -------
    traj : dict
    TODO: traj description
    TODO: remove `ts`?
    """

    traj = {}
    natoms = getNatoms(lines)
    nlines = natoms + header
    firstline = ((tsnum-1)*nlines)
    lastline = firstline + nlines-1
    for line in range(firstline, lastline):
        if lines[line].startswith("ITEM: TIMESTEP"):
            ts = int(lines[line+1])
            traj[ts] = {}
            traj[ts]["atom"] = {}
        if any(s in lines[line] for s in ["BOX BOUNDS", "xy xz yz pp pp"]):
            boxx = list(map(float, lines[line+1].split()))
            boxy = list(map(float, lines[line+2].split()))
            boxz = list(map(float, lines[line+3].split()))
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            traj[ts]["boxx"] = (boxx)
            traj[ts]["boxy"] = (boxy)
            traj[ts]["boxz"] = (boxz)
            traj[ts]["boxsize"] = (lx, ly, lz)
        if lines[line].startswith("ITEM: ATOMS"):
            labels = lines[line].split()
            for atomid in range(1, natoms + 1):
                atomvals = lines[line+atomid].split()
                traj[ts]["atom"][atomid] = {}
                for label_idx, k in enumerate(labels[3:], start=1):
                    if k == 'type':
                        traj[ts]["atom"][atomid][k] = int(atomvals[label_idx])
                    else:
                        traj[ts]["atom"][atomid][k] = float(
                            atomvals[label_idx])
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
        A trajectory dict with more than one trajectory. Keys are the
        timestep labels.
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


def getAtomData(filename):
    """
    Reads in a data file, returns atom data, including mass.

    Parameters
    ----------
    filename : str
        Name of the datafile to be read.

    Returns
    -------
    atomdata : dict
        Each entry of the dictionary takes the form:
        atomid (int): {'mol': int,
                       'type': str,
                       'charge': float,
                       'x': float,
                       'y': float,
                       'z': float}

    masses : dict
        The atomic masses. Each entry takes the form:
        type (str): mass (float)
    """

    lines = readAll(filename)
    atomnames = getAtomType(filename)
    atomdata = {}
    masses = {}
    for line_idx, x in enumerate(lines):
        if "atoms" in lines[line_idx]:
            natoms = int(lines[line_idx].split("atoms")[0])

        if "atom types" in lines[line_idx]:
            ntypes = int(lines[line_idx].split("atom")[0])
            print("ntypes ", ntypes)

        if "Masses" in lines[line_idx]:
            for atomtype in range(int(ntypes)):
                el = lines[line_idx + 2 + atomtype].split()
                print(el)
                masses[atomnames[int(el[0])]] = float(el[1])

        if "Atoms" in lines[line_idx]:
            columns = ['mol', 'type', 'charge', 'x', 'y', 'z']
            print(columns)
            print(columns[1])
            for atom in range(int(natoms)):
                el = lines[line_idx + 2 + atom].split()
                atomid = int(el[0])
                atomdata[atomid] = {}
                for column in columns:
                    if column in ['mol', 'type']:
                        atomdata[atomid][column] = \
                            int(el[columns.index(column) + 1])
                    else:
                        atomdata[atomid][column] = \
                            float(el[columns.index(column) + 1])
                    # print(atomnames)
                atomdata[atomid]['type'] = \
                    atomnames[int(atomdata[atomid]['type'])]
    return atomdata, masses
