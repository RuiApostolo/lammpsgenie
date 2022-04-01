import re
import lmptools.readfiles as rdfl


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
            print()
            print(f"Found {natomtypes} atom types.")
            print()
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


def getAtomData(filename):
    """
    Reads in a data file, returns reduced atom data, and masses.

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

    lines = rdfl.readAll(filename)
    atomnames = getAtomType(filename)

    for line_idx, line in enumerate(lines):
        if "atoms" in line:
            natoms = _getFirst(line, "atoms")

        elif "atom types" in line:
            natomtypes = _getFirst(line, "atom types")
            print("ntypes ", natomtypes)

        elif "Masses" in line:
            masses = _getMasses(lines, line_idx, natomtypes, atomnames)

        elif "Atoms" in line:
            atomdata = _getAtoms(lines, line_idx, natoms, atomnames)

    return atomdata, masses


# TODO
def getAllAtomData(filename):
    lines = rdfl.readAll(filename)
    atomnames = getAtomType(filename)

    bonddata = {}
    angledata = {}
    dhdata = {}
    impdata = {}
    boxsize = {}
    pair = {}
    bond = {}
    angle = {}
    dh = {}
    imp = {}

    numdata = {}
    limits = ['xlo xhi', 'ylo yhi', 'zlo zhi']

    for line_idx, line in enumerate(lines):
        if "atoms" in line:
            natoms = _getFirst(line, "atoms")
            numdata['natoms'] = natoms
        elif "bonds" in line:
            nbonds = _getFirst(line, "bonds")
            numdata['nbonds'] = nbonds
        elif "angles" in line:
            nangles = _getFirst(line, "angles")
            numdata['nangles'] = nangles
        elif "dihedrals" in line:
            ndh = _getFirst(line, "dihedrals")
            numdata['ndh'] = ndh
        elif "impropers" in line:
            nimp = _getFirst(line, "impropers")
            numdata['nimp'] = nimp
        elif "atom types" in line:
            natomtypes = _getFirst(line, "atom types")
            numdata['natomtypes'] = natomtypes
            print("natomtypes "+str(natomtypes))
        elif "bond types" in line:
            nbondtypes = _getFirst(line, "bond")
            numdata['nbondtypes'] = nbondtypes
            print("nbondtypes "+str(nbondtypes))
        elif "angle types" in line:
            nangletypes = _getFirst(line, "angle")
            numdata['nangletypes'] = nangletypes
            print("nangletypes "+str(nangletypes))
        elif "dihedral types" in line:
            ndhtypes = _getFirst(line, "dihedral")
            numdata['ndhtypes'] = ndhtypes
            print("ndhtypes "+str(ndhtypes))
        elif "improper types" in line:
            nimptypes = _getFirst(line, "improper")
            numdata['nimptypes'] = nimptypes
            print("nimptypes "+str(nimptypes))

        elif any(limit in line for limit in limits):
            low, high = line.split()[2:4]
            boxsize[low], boxsize[high] = map(float, line.split()[0:2])

        elif "Pair Coeffs" in line:
            columns = ['eps', 'sig', '#', 'label']
            for j in range(int(natomtypes)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        paircoeffid = int(el[0])
                        pair[paircoeffid] = {}
                        pair[paircoeffid][str(col)] = {}
                    pair[paircoeffid][str(col)] = (el[col_idx + 1])
            print("Pair Coeffs:")
            print(pair)

        elif "Bond Coeffs" in line:
            columns = ['k', 'r', '#', 'label']
            for j in range(int(nbondtypes)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        bondcoeffid = int(el[0])
                        bond[bondcoeffid] = {}
                        bond[bondcoeffid][str(col)] = {}
                    bond[bondcoeffid][str(col)] = (el[col_idx+1])
            print("Bond Coeffs")
            print(bond)

        elif "Angle Coeffs" in line:
            columns = ['k', 'theta', '#', 'label']
            for j in range(int(nangletypes)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        anglecoeffid = int(el[0])
                        angle[anglecoeffid] = {}
                        angle[anglecoeffid][str(col)] = {}
                    angle[anglecoeffid][str(col)] = (el[col_idx+1])
            print("Angle Coeffs")
            print(angle)

        elif "Dihedral Coeffs" in line:
            columns = ['a', 'b', 'c', 'd', '#', 'label']
            for j in range(int(ndhtypes)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        dhcoeffid = int(el[0])
                        dh[dhcoeffid] = {}
                        dh[dhcoeffid][str(col)] = {}
                    dh[dhcoeffid][str(col)] = (el[col_idx+1])
            print("Dihedral Coeffs")
            print(dh)

        elif "Improper Coeffs" in line:
            columns = ['a', 'b', 'c', 'd', '#', 'label']
            for j in range(int(nimptypes)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        impcoeffid = int(el[0])
                        imp[impcoeffid] = {}
                        imp[impcoeffid][str(col)] = {}
                    imp[impcoeffid][str(col)] = (el[col_idx+1])
            print("Improper Coeffs")
            print(imp)

        elif "Masses" in line:
            masses = _getMasses(lines, line_idx, natomtypes, atomnames)

        elif "Atoms" in line:
            atomdata = _getAtoms(lines, line_idx, natoms, atomnames)

        elif "Bonds" in line:
            columns = ['bondtype', 'at1', 'at2']
            for j in range(int(nbonds)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        bondid = int(el[0])
                        bonddata[bondid] = {}
                        bonddata[bondid][str(col)] = {}
                    bonddata[bondid][str(col)] = int(el[col_idx+1])

        elif "Angles" in line:
            columns = ['angletype', 'at1', 'at2', 'at3']
            for j in range(int(nangles)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        angleid = int(el[0])
                        angledata[angleid] = {}
                        angledata[angleid][str(col)] = {}
                    angledata[angleid][str(col)] = int(el[col_idx+1])

        elif "Dihedrals" in line:
            columns = ['dhtype', 'at1', 'at2', 'at3', 'at4']
            for j in range(int(ndh)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        dhid = int(el[0])
                        dhdata[dhid] = {}
                        dhdata[dhid][str(col)] = {}
                    dhdata[dhid][str(col)] = int(el[col_idx+1])

        elif "Impropers" in line:
            columns = ['imptype', 'at1', 'at2', 'at3', 'at4']
            for j in range(int(nimp)):
                el = lines[line_idx+2+j].split()
                for col_idx, col in enumerate(columns):
                    if col_idx == 0:
                        impid = int(el[0])
                        impdata[impid] = {}
                        impdata[impid][str(col)] = {}
                    impdata[impid][str(col)] = int(el[col_idx+1])

    return atomdata, bonddata, angledata, dhdata, impdata, masses, boxsize, \
        numdata, pair, bond, angle, dh, imp


###############################################################################
#                             Protected Functions                             #
###############################################################################
def _getFirst(line, word):
    """
    Returns split part before 'word'

    Parameters
    ----------
    line : str
        Line to be split
    word : str
        Word to split line on

    Returns
    -------
    int
        Returns integer before word.
    """

    return int(line.split(word)[0])


def _getMasses(lines, line_idx, natomtypes, atomnames):
    """
    Returns split part before 'word'

    Parameters
    ----------
    lines : list of str
        List of lines from dumpfile.
    line_idx : int
        Index of current line.
    natomtypes : int
        Number of atom types.
    atomnames : dict
        A dictionary with the LAMMPS atom type numbers as keys, and the
        custom atom type names as values. From getAtomType().

    Returns
    -------
    masses : dict
        The atomic masses. Each entry takes the form:
        type (str): mass (float)
    """

    masses = {}
    for atomtype in range(int(natomtypes)):
        el = lines[line_idx + 2 + atomtype].split()
        print(el)
        # assign masses by atom type
        masses[atomnames[int(el[0])]] = float(el[1])
    return masses


def _getAtoms(lines, line_idx, natoms, atomnames):
    """
    Reads 'Atoms' data section in a LAMMPS data file, and returns a dict
    of the data, with `atomid` as keys.

    Parameters
    ----------
    lines : list of str
        List of lines from dumpfile.
    line_idx : int
        Index of current line.
    natoms : int
        Number of atoms.
    atomnames : dict
        A dictionary with the LAMMPS atom type numbers as keys, and the
        custom atom type names as values. From getAtomType().

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
    """

    atomdata = {}
    columns = ['mol', 'type', 'charge', 'x', 'y', 'z']
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
    return atomdata
