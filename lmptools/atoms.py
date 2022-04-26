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
            #  print("ntypes ", natomtypes)

        elif "Masses" in line:
            masses = _getMasses(lines, line_idx, natomtypes, atomnames)

        elif "Atoms" in line:
            atomdata = _getAtoms(lines, line_idx, natoms, atomnames)

    return atomdata, masses


# TODO docstring
def getAllAtomData(filename): #noqa C901
    lines = rdfl.readAll(filename)
    atomnames = getAtomType(filename)

    boxsize = {}

    numdata = {}
    limits = ['xlo xhi', 'ylo yhi', 'zlo zhi']
    properties = {
        'atoms': 0,
        'bonds': 0,
        'angles': 0,
        'dihedrals': 0,
        'impropers': 0,
        'atom types': 0,
        'bond types': 0,
        'angle types': 0,
        'dihedral types': 0,
        'improper types': 0
    }
    coeffs = {'Pair Coeffs': ['atom types', {}],
              'Bond Coeffs': ['bond types', {}],
              'Angle Coeffs': ['angle types', {}],
              'Dihedral Coeffs': ['dihedral types', {}],
              'Improper Coeffs': ['improper types', {}]
              }
    badis = {'Bonds': {},
             'Angles': {},
             'Dihedrals': {},
             'Impropers': {}
             }

    for line_idx, line in enumerate(lines):
        if any(propert in line for propert in properties):
            for propert in properties:
                if propert in line:
                    properties[propert] = _getFirst(line, propert)

        elif any(limit in line for limit in limits):
            low, high = line.split()[2:4]
            boxsize[low], boxsize[high] = map(float, line.split()[0:2])

        elif any(coeff in line for coeff in coeffs):
            for coeff in coeffs:
                if coeff in line:
                    coeffs[coeff][1] = _getCoeffs(lines,
                                                  line_idx,
                                                  properties[coeffs[coeff][0]])

        elif "Masses" in line:
            masses = _getMasses(lines,
                                line_idx,
                                properties['atom types'],
                                atomnames)

        elif "Atoms" in line:
            atomdata = _getAtoms(lines,
                                 line_idx,
                                 properties['atoms'],
                                 atomnames)

        elif any(badi in line for badi in badis):
            for badi in badis:
                if badi in line:
                    badis[badi] = _getBADI(lines,
                                           line_idx,
                                           properties[badi.lower()])

    pair, bond, angle, dihedral, imp = [coeffs[coeff][1] for coeff in coeffs]
    bonddata, angledata, dhdata, impdata = [badis[badi] for badi in badis]

    return atomdata, bonddata, angledata, dhdata, impdata, masses, boxsize, \
        numdata, pair, bond, angle, dihedral, imp


def getAtomRange(atomdata, atom1, atom2):
    """
    Returns list of atom ids between atom1 and atom2, inclusive.

    Parameters
    ----------
    atomdata : dict
        Each entry of the dictionary takes the form:
        atomid (int): {'mol': int,
                       'type': str,
                       'charge': float,
                       'x': float,
                       'y': float,
                       'z': float}
    atom1 : int
        Index of first atom to select.
    atom2 : int
        Index of last atom to select.

    Returns
    -------
    atomrange : list of ints
        Returns list of atom ids.
    """
    sortedatomnum = sorted(atomdata.keys())
    try:
        atom1 = 0 if atom1 == 'first' else sortedatomnum.index(atom1)
        atom2 = None if atom2 == 'last' else sortedatomnum.index(atom2)
        if atom2 is None:
            atrange = sortedatomnum[atom1:None]
        else:
            atrange = sortedatomnum[atom1:atom2 + 1]
        return atrange
    except(ValueError):
        return []


def getAtomsByType(atomdata, *types):
    """
    Returns list of atom ids that correspond to `types`.

    Parameters
    ----------
    atomdata : dict

        Each entry of the dictionary takes the form:
        atomid (int): {'mol': int,
                       'type': str,
                       'charge': float,
                       'x': float,
                       'y': float,
                       'z': float}
    *types: str
        Atom types in string form.

    Returns
    -------
    atomrange : list of ints
        Returns list of atom ids.
    """

    atomrange = []
    sortedatomnum = sorted(atomdata.keys())
    for atom in sortedatomnum:
        for typ in types:
            if atomdata[atom]['type'] == typ:
                atomrange.append(atom)
    return atomrange


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
        # replace type from int to str
        atomdata[atomid]['type'] = \
            atomnames[int(atomdata[atomid]['type'])]
    return atomdata


def _getParams(params, ptype):
    """
    Returns split parameters as floats, not including index.

    Parameters
    ----------
    params : list of str
        Parameters to be split.

    Returns
    -------
    params : list of floats
        a list of the parameters transformed to floats.
    """

    params = params.split()
    params = list(map(ptype, params))
    return params[1:]


def _getCoeffs(lines, line_idx, number):
    """
    Returns dictionary with coefficients.

    Parameters
    ----------
    lines : list of str
        List of lines from dumpfile.
    line_idx : int
        Index of current line
    number : int
        Number of lines to read through - natomtypes, nbondtypes, etc.

    Returns
    -------
    dict
        Dictionary where keys are the type number (int), and the values are
        lists with the parameters of the pair, bond, angle, dihedral, or
        improper.
    """

    result = {}
    for coeff in range(1, int(number) + 1):
        params, comment = lines[line_idx + 1 + coeff].split('#')
        result[coeff] = _getParams(params, float) + [comment]
    #  print("Angle Coeffs")
    #  print(angle)
    return result


# TODO: get better name
def _getBADI(lines, line_idx, number):
    """
    Returns dictionary with atom numbers for the bond/angle/dihedral/improper.

    Parameters
    ----------
    lines : list of str
        List of lines from dumpfile.
    line_idx : int
        Index of current line
    number : int
        Number of lines to read through - nbonds, nangles, etc.

    Returns
    -------
    dict
        Dictionary where keys are the index (int), and the values are lists
        of the atoms in the bond, angle, dihedral, or improper.
    """
    result = {}
    for line in range(1, int(number) + 1):
        result[line] = _getParams(lines[line_idx + 1 + line], int)
    return result
