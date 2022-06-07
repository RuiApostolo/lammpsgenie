"""
Function that gets atomic information from already loaded data and dump files.
"""

import lammpsgenie.readfiles as rdfl
from re import compile, match
from collections import defaultdict


__all__ = [
    'getNatoms',
    'getAtomType',
    'getAtomData',
    'getAllAtomData',
    'getAtomRange',
    'getAtomsByType',
    'getTotalMass',
    'getCOM',
]


g_coeffs = {'Pair Coeffs': 'atom types',
            'Bond Coeffs': 'bond types',
            'Angle Coeffs': 'angle types',
            'Dihedral Coeffs': 'dihedral types',
            'Improper Coeffs': 'improper types',
            }


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

    for lineindex, _ in enumerate(lines):
        if lines[lineindex].startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[lineindex+1])
            return natoms
    return None


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

    p = compile(r"(\s*\d*\s*)(?P<name>\w+)(\s*\w*)")
    with open(filename, "r") as datafile:
        lines = datafile.read().splitlines()
    for lineindex, line in enumerate(lines):
        if "atom types" in line:
            natomtypes = int(line.split()[0])
            atomnames = {}
        if "Pair Coeffs" in line:
            for i in range(1, natomtypes + 1):
                matches = match(p, lines[lineindex + i + 1].split("#")[1])
                try:
                    atomnames[i] = matches.group('name')
                except(AttributeError):
                    atomnames[i] = str(i)
            return atomnames
    return {}


def getAtomData(filename):
    """
    Retrives reduced atom data, and masses, from a LAMMPS data file.

    Parameters
    ----------
    filename : str
        Name of the datafile to be read.

    Returns
    -------
    Topology : dict, comprised of
        atomnames : dict
            A dictionary with the LAMMPS atom type numbers as keys, and the
            custom atom type names as values. From getAtomType().

        masses : dict
            The atomic masses. Each entry takes the form:
            type (str): mass (float)

        atomdata : dict
            Each entry of the dictionary takes the form:
            atomid (int): {'mol': int,
                           'type': int,
                           'charge': float,
                           'x': float,
                           'y': float,
                           'z': float}
    """

    lines = rdfl.readAll(filename)
    atomnames = getAtomType(filename)

    for line_idx, line in enumerate(lines):
        if "atoms" in line:
            natoms = _getFirst(line, "atoms")

        elif "atom types" in line:
            natomtypes = _getFirst(line, "atom types")

        elif "Masses" in line:
            masses = _getMasses(lines, line_idx, natomtypes, atomnames)

        elif "Atoms" in line:
            atomdata = _getAtoms(lines, line_idx, natoms, atomnames)

    return {
        'topologycounts': {},
        'boxsize': {},
        'atomnames': atomnames,
        'pairtypes': {},
        'bondtypes': {},
        'angletypes': {},
        'dihedraltypes': {},
        'impropertypes': {},
        'masses': masses,
        'atomdata': atomdata,
        'bonddata': {},
        'angledata': {},
        'dihedraldata': {},
        'improperdata': {},
            }


def getAllAtomData(filename):  # noqa C901
    """
    Retrieves every piece of atom data from a LAMMPS data file.

    Parameters
    ----------
    filename : str
        Name of the datafile to be read.

    Returns
    -------
    topology : dict, comprised of
        topologycounts : dict
            A dictionary that contains the number, and number of types of
            atoms, bonds, angles, dihedrals, and impropers.
            Each entry takes the form:
            property (str) : int
            where property can be one of:
            'atoms', 'bonds', 'angles', 'dihedrals', 'impropers', 'atom types',
            'bond types', 'angle types', 'dihedral types' and 'improper types'.

        atomnames : dict
            A dictionary with the LAMMPS atom type numbers as keys, and the
            custom atom type names as values. From getAtomType().

        boxsize : dict
            The simulaton box size. The dictionary takes the form:
            coordinate (str) : float
            where coordinate is can take the values:
            'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi'
            in accordance with LAMMPS data file syntax.

        pairtypes, bondtypes, angletypes, dihedraltypes, impropertypes : dict
            The pair, bond, angle, dihedral, and improper coefficients.
            Each dictionary entry takes the form:
            property_type (int) : [*coeffs, 'comment'] list of floats + str
            The number of coefficients can vary depending on the particular
            property. Pair, bond, angles, and impropers usually have 2,
            dihedrals usually have 4. The function accepts any number.

        masses : dict
            The atomic masses. Each entry takes the form:
            type (str): mass (float)

        atomdata : dict
            Each entry of the dictionary takes the form:
            atomid (int): {'mol': int,
                           'type': int,
                           'charge': float,
                           'x': float,
                           'y': float,
                           'z': float}

        bonddata : dict
            The bond data. Each entry of the dictionary takes the form:
            bondid (int) : [bondtype, atom1, atom2] (list of ints)

        angledata : dict
            The angle data. Each entry of the dictionary takes the form:
            angleid (int) : [angletype, atom1, atom2, atom3] (list of ints)

        dihedraldata : dict
            The dihedral data. Each entry of the dictionary takes the form:
            dihedralid (int) : [dihedraltype, atom1, atom2, atom3, atom4]
                               (list of ints)

        impdata : dict
            The improper data. Each entry of the dictionary takes the form:
            improperid (int) : [impropertype, atom1, atom2, atom3]
                               (list of ints)
    """

    lines = rdfl.readAll(filename)
    atomnames = getAtomType(filename)

    boxsize = {}

    limits = ['xlo xhi', 'ylo yhi', 'zlo zhi']
    topologycounts = {
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
    coeffs = {coeff: {'name': value, 'list': {}}
              for coeff, value in g_coeffs.items()}
    badis = {'Bonds': {},
             'Angles': {},
             'Dihedrals': {},
             'Impropers': {}
             }

    for line_idx, line in enumerate(lines):
        if any(topologycount in line for topologycount in topologycounts):
            for topologycount in topologycounts:
                if topologycount in line:
                    topologycounts[topologycount] = \
                        _getFirst(line, topologycount)

        elif any(limit in line for limit in limits):
            low, high = line.split()[2:4]
            boxsize[low], boxsize[high] = map(float, line.split()[0:2])

        elif any(coeff in line for coeff in coeffs):
            for coeff in coeffs:
                if coeff in line:
                    coeffs[coeff]['list'] = \
                        _getCoeffs(lines,
                                   line_idx,
                                   topologycounts[coeffs[coeff]['name']])

        elif "Masses" in line:
            masses = _getMasses(lines,
                                line_idx,
                                topologycounts['atom types'],
                                atomnames)

        elif "Atoms" in line:
            atomdata = _getAtoms(lines,
                                 line_idx,
                                 topologycounts['atoms'],
                                 atomnames)

        elif any(badi in line for badi in badis):
            for badi in badis:
                if badi in line:
                    badis[badi] = _getBADI(lines,
                                           line_idx,
                                           topologycounts[badi.lower()])

    pairtypes, bondtypes, angletypes, dihedraltypes, impropertypes = \
        [coeffs[coeff]['list'] for coeff in coeffs]
    bonddata, angledata, dihedraldata, improperdata = \
        [badis[badi] for badi in badis]

    return {
        'topologycounts': topologycounts,
        'boxsize': boxsize,
        'atomnames': atomnames,
        'pairtypes': pairtypes,
        'bondtypes': bondtypes,
        'angletypes': angletypes,
        'dihedraltypes': dihedraltypes,
        'impropertypes': impropertypes,
        'masses': masses,
        'atomdata': atomdata,
        'bonddata': bonddata,
        'angledata': angledata,
        'dihedraldata': dihedraldata,
        'improperdata': improperdata,
            }


def getAtomRange(atomdata, atom1, atom2):
    """
    Returns list of atom ids between atom1 and atom2, inclusive.

    Parameters
    ----------
    atomdata : dict
        Each entry of the dictionary takes the form:
        atomid (int): {'mol': int,
                       'type': int,
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


def getAtomsByType(atomdata, atomnames, *types):
    """
    Returns list of atom ids that correspond to `types`.

    Parameters
    ----------
    atomdata : dict
        Each entry of the dictionary takes the form:
        atomid (int): {'mol': int,
                       'type': int,
                       'charge': float,
                       'x': float,
                       'y': float,
                       'z': float}

    atomnames : dict
        A dictionary with the LAMMPS atom type numbers as keys, and the
        custom atom type names as values. From getAtomType().

    types : str
        Atom types in string form, like 'C', 'H', or 'SCP'.

    Returns
    -------
    atomrange : list of ints
        Returns list of atom ids.
    """

    atomrange = []
    atomnames = _invertDictToList(atomnames)
    sortedatomnum = sorted(atomdata.keys())
    for atom in sortedatomnum:
        for typ in types:

            if atomdata[atom]['type'] in [a for a in atomnames[typ]]:
                atomrange.append(atom)
    return atomrange


def getTotalMass(traj, masses):
    """
    Returns total mass of atoms in first timestep of the provided
    trajectory.

    Parameters
    ----------
    traj : dict
        A trajectory dict with more than one timestep. Keys are the
        timestep labels.

    masses : dict
        The atomic masses. Each entry takes the form:
        type (str): mass (float)

    Returns
    -------
    totalmass : float
        The sum of masses of all atoms in the timestep.
    """

    frame = sorted(traj.keys())[0]
    totalmass = 0.0
    # for key in dict:
    for atom in traj[frame]["atom"]:
        totalmass += masses[traj[frame]["atom"][atom]['type']]
    return totalmass


def getCOM(traj, tsrange, atomlist, masses):
    """
    Returns centre of mass position for chosen atoms in timestep range.

    Parameters
    ----------
    traj : dict
        A trajectory dict with more than one timestep. Keys are the
        timestep labels.
    tsrange : list of int
        List of timestep numbers.
    atomlist : list of ints
        List of atom ids.
    masses : dict
        The atomic masses. Each entry takes the form:
        type (str): mass (float)

    Returns
    -------
    com : dict
        timestep label (int) : (x, y, z) (tuple of floats)
    """

    com = {}
    totalmass = getTotalMass(traj, masses)
    for ts in tsrange:
        xt, yt, zt = 0.0, 0.0, 0.0
        for atom in [a for a in traj[ts]["atom"].keys() if a in atomlist]:
            mass = masses[traj[ts]["atom"][atom]['type']]
            xt += traj[ts]["atom"][atom]['x'] * mass
            yt += traj[ts]["atom"][atom]['y'] * mass
            zt += traj[ts]["atom"][atom]['z'] * mass
        com[ts] = (xt/totalmass, yt/totalmass, zt/totalmass)
    return com


###############################################################################
#                             Protected Functions                             #
###############################################################################
def _getFirst(line, word):
    """
    Returns an integer in 'line' before a specified delimiter 'word'.

    Parameters
    ----------
    line : str
        Line to be split.
    word : str
        Delimiter used to split 'line'.

    Returns
    -------
    int
        Returns integer before 'word'.
    """

    return int(line.split(word)[0])


def _getMasses(lines, line_idx, natomtypes, atomnames):
    """
    Returns masses for all atoms in lines.

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
        word = lines[line_idx + 2 + atomtype].split()
        # assign masses by atom type
        masses[atomnames[int(word[0])]] = float(word[1])
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
                       'type': int,
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
        # replace type from int to str -- BAD IDEA, can't then demux back to
        # ints if more than one merged file shares same atom type names.
        #  atomdata[atomid]['type'] = \
        #      atomnames[int(atomdata[atomid]['type'])]
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
    return result


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


def _invertDictToList(dictionary):
    # from a dict of {key (int): value (str)}, inverts to:
    # {value (str): [key, key, key, key] (list of ints)}
    # each new_key will have new_values that had common old_values
    reverse_dict = defaultdict(list)
    for key, value in dictionary.items():
        reverse_dict[value].append(key)
    return reverse_dict
