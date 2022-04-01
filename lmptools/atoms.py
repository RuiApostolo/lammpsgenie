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

    lines = rdfl.readAll(filename)
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
