#!/bin/env python3
<<<<<<< HEAD:lmptools/mergedatafiles.py
"""
Functions and script required to merge LAMMPS data diles.
"""

import lammpsgenie.atoms as atoms
from lammpsgenie._version import __version__
from lammpsgenie._name import _name
>>>>>>> master:lammpsgenie/mergedatafiles.py
from yaml import full_load
from sys import argv
from copy import deepcopy
from datetime import datetime
from os.path import basename
from warnings import warn
import operator


__all__ = [
    'readInputFile',
    'readTopology',
    'readTopologies',
    'limitsTopology',
    'limitsAllTopologies',
    'absoluteLimitsTopologies',
    'shiftTopology',
    'shiftTopologies',
    'mergeTopologies',
    'writeTopology',
]

coeffs_equiv = {
    'Pair Coeffs': 'pairtypes',
    'Bond Coeffs': 'bondtypes',
    'Angle Coeffs': 'angletypes',
    'Dihedral Coeffs': 'dihedraltypes',
    'Improper Coeffs': 'impropertypes',
    }


class MissingSettingsFile(IOError):
    pass


class TooManyArguments(IOError):
    pass


class ValueExists(ValueError):
    pass


args = argv
_this_file = basename(__file__)


def readInputFile(args):
    """
    Reads settings file.

    Parameters
    ----------
    args : list of str
        Command-line arguments passed. Accepts name of YAML settings file
        or None. If None, it will search the execution directory for
        'merge.yaml' and 'merge.yml' and use the first one it finds.

    Returns
    -------
    settings : dict of dicts
        Each dict should contain:
        'filename' :
            'minormax': function
            'value': float

    Raises
    ------
    MissingSettingFile
        If passed no settings file name, and the defaults aren't found.

    TooManyArguments
        If passed more arguments than required.
    """

    # trims name of file run (argv[0])
    args = args[1:]
    try:
        if len(args) < 1:
            for filename in ['merge.yaml', 'merge.yml']:
                try:
                    return _openYaml(filename)

                except IOError:
                    continue
            else:
                raise MissingSettingsFile
        elif len(args) > 1:
            raise TooManyArguments
        else:
            filename = args[0]
            return _openYaml(filename)

    except MissingSettingsFile as exception:
        message = "Missing settings file: merge.yaml or merge.yml"
        print(message)
        raise MissingSettingsFile(message) from exception

    except TooManyArguments as exception:
        message = "This script takes only one argument, the settings file."
        print(message)
        raise MissingSettingsFile(message) from exception


def readTopology(file):
    """
    Reads Topology from settings dict.

    Parameters
    ----------
    file : dict
        Name of the datafile to be read.

    Returns
    -------
    topology : dict
        Dictionary in the style of atoms.getAllAtomData()
    """

    return atoms.getAllAtomData(file)


def readTopologies(settings):
    """
    Reads Topologies from settings dict.

    Parameters
    ----------
    settings : list of dicts
        Settings, from readInputFile().

    Returns
    -------
    topologies : nested dicts
        Each top level entry has the shape:
        'filename': Dictionary in the style of atoms.getAllAtomData()
    """

    topologies = {}
    for filename in settings:
        topologies[filename] = readTopology(filename)
    return topologies


def limitsTopology(topology):
    """
    Finds the min and max values of the coordinate of a topology file.

    Parameters
    ----------
    topology : nested dicts
        dictionary of the form:
        'filename': Dictionary in the style of atoms.getAllAtomData()

    Returns
    -------
    limits : dict
        Takes the form:
        {min: {'x': float, 'y': float, 'z': float},
         max: {'x': float, 'y': float, 'z': float}}

    Notes
    -----
    Uses the actual functions min() and max() and not str.
    """

    limits = {min: {}, max: {}}
    coords = ['x', 'y', 'z']
    for limit in limits:
        for coord in coords:
            limits[limit][coord] = _limit(limit, coord, topology)
    return limits


def limitsAllTopologies(topologies):
    """
    Finds the min and max values of each coordinate for each topology file.

    Parameters
    ----------
    topologies : nested dicts
        Each top level entry has the shape:
        'filename': Dictionary in the style of atoms.getAllAtomData()

    Returns
    -------
    limits : nested dict
        Each entry takes the form:
        filename: {
            min: {'x': float, 'y': float, 'z': float},
            max: {'x': float, 'y': float, 'z': float}
                   }

    Notes
    -----
    Uses the actual functions min() and max() and not str.
    """

    all_limits = {}
    for topology in topologies:
        all_limits[topology] = limitsTopology(topologies[topology])
    return all_limits


def absoluteLimitsTopologies(topologies):
    """
    Finds the min and max values of each coordinate for a series of topologies.
    Does not save which topology each value originated from.

    Parameters
    ----------
    topologies : nested dicts
        Each top level entry has the shape:
        'filename': Dictionary in the style of atoms.getAllAtomData()

    Returns
    -------
    limits : dict
        Takes the form:
        {min: {'x': float, 'y': float, 'z': float},
         max: {'x': float, 'y': float, 'z': float}}

    Notes
    -----
    Uses the actual functions min() and max() and not str.
    """

    limits = {min: {'x': float("+inf"),
                    'y': float("+inf"),
                    'z': float("+inf")},
              max: {'x': float("-inf"),
                    'y': float("-inf"),
                    'z': float("-inf")}}
    for topology in topologies:
        this_topology_limits = limitsTopology(topologies[topology])
        for limit in this_topology_limits:
            for coord in this_topology_limits[limit]:
                limits[limit][coord] = \
                    limit(limits[limit][coord],
                          this_topology_limits[limit][coord])
    return limits


def shiftTopology(topology, limits, settings, axis='z'):
    """
    Moves the coordinates of the topology file to the min/max along a given
    axis, as chosen in the settings file.

    Parameters
    ----------
    topology : nested dicts
        dictionary of the form:
        'filename': Dictionary in the style of atoms.getAllAtomData()

    limits : dict
        Takes the form:
        {min: {'x': float, 'y': float, 'z': float},
         max: {'x': float, 'y': float, 'z': float}}

    settings : dict
        Each dict should contain:
        'filename' :
            'minormax': function
            'value': float

    axis : {'z', 'x', 'y'}
        The axis along which to apply the coordinate shift.

    Returns
    -------
    topology : nested dicts
        Dictionary in the style of atoms.getAllAtomData(), with shifted
        coordinates.
    """

    # copy dictionary. d1 = d2.copy() creates only a pointer
    new_topology = deepcopy(topology)
    # what to add to current coordinates
    # 'value' - min/max(minimum_coord, maximum_coord)
    delta = settings['value'] - \
        settings['minormax'](limits[min][axis], limits[max][axis])
    for atom in new_topology['atomdata']:
        new_topology['atomdata'][atom][axis] += delta
    return new_topology


def shiftTopologies(topologies, limits, settings, axis='z'):
    """
    Moves the coordinates of the topology file to the min/max along a given
    axis, as chosen in the settings file.

    Parameters
    ----------
    topologies : nested dicts
        Each top level entry has the shape:
        'filename': Dictionary in the style of atoms.getAllAtomData()

    limits : nested dict
        Each entry takes the form:
        filename: {
            min: {'x': float, 'y': float, 'z': float},
            max: {'x': float, 'y': float, 'z': float}
                   }

    settings : dict of dicts
        Each dict should contain:
        'filename' :
            'minormax': function
            'value': float

    axis : {'z', 'x', 'y'}
        The axis along which to apply the coordinate shift.

    Returns
    -------
    topologies : nested dicts
        Each top level entry has the shape:
        'filename': Dictionary in the style of atoms.getAllAtomData()
    """

    new_topologies = {}
    for topology in topologies:
        new_topologies[topology] = shiftTopology(
                                    topologies[topology],
                                    limits[topology],
                                    settings[topology],
                                    axis)
    return new_topologies


def mergeTopologies(topologies, newboxsize):
    """
    Merges a dictionary containing several topologies into one single
    topology.

    Parameters
    ----------
    topologies : nested dicts
        Each top level entry has the shape:
        'filename': Dictionary in the style of atoms.getAllAtomData()

    newboxsize : dict
        New simulation box size, intended to be read from settings file.
        The dictionary takes the form:
        coordinate (str) : float
        where `coordinate` can have the values:
        'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi'
        in accordance with LAMMPS data file syntax.

    Returns
    -------
    topology : nested dicts
        New, merged topology. Dictionary in the style of
        atoms.getAllAtomData().
    """

    merged_topology = {}
    property_pairs = {
        'atomnames':
            {'key': 'atom types',
             'listid': None,
             'listelements': None},
        'pairtypes':
            {'key': 'atom types',
             'listid': None,
             'listelements': None},
        'bondtypes':
            {'key': 'bond types',
             'listid': None,
             'listelements': None},
        'angletypes':
            {'key': 'angle types',
             'listid': None,
             'listelements': None},
        'dihedraltypes':
            {'key': 'dihedral types',
             'listid': None,
             'listelements': None},
        'impropertypes':
            {'key': 'improper types',
             'listid': None,
             'listelements': None},
        'bonddata':
            {'key': 'bonds',
             'listid': 'bond types',
             'listelements': 'atoms'},
        'angledata':
            {'key': 'angles',
             'listid': 'angle types',
             'listelements': 'atoms'},
        'dihedraldata':
            {'key': 'dihedrals',
             'listid': 'dihedral types',
             'listelements': 'atoms'},
        'improperdata':
            {'key': 'impropers',
             'listid': 'improper types',
             'listelements': 'atoms'},
    }
    molid = 1
    for topology in topologies:
        # no changes on first topology
        if len(merged_topology) == 0:
            merged_topology = deepcopy(topologies[topology])
            # reset molid just in case
            merged_topology['atomdata'] = _setMolid(
                merged_topology['atomdata'], molid)
        else:
            shift_values = merged_topology['topologycounts']
            shift_values[None] = None
            for propert in property_pairs:
                merged_topology[propert] = _mergeDicts(
                    merged_topology[propert], _shiftDictOfLists(
                        topologies[topology][propert],
                        shift_values[property_pairs[propert]['key']],
                        shift_values[property_pairs[propert]['listid']],
                        shift_values[property_pairs[propert]['listelements']],
                        molid))
            # atomdata
            merged_topology['atomdata'] = _mergeDicts(
                merged_topology['atomdata'], _shiftDictOfDicts(
                    topologies[topology]['atomdata'],
                    shift_values['atoms'],
                    'type',
                    shift_values['atom types'],
                    molid))
            # masses
            merged_topology['masses'] = _mergeDicts(
                merged_topology['masses'], topologies[topology]['masses'])
            # add topology counts
            merged_topology['topologycounts'] = _combineDicts(
                shift_values, topologies[topology]['topologycounts'])
        molid += 1
    # limits/boxsize
    merged_topology['boxsize'] = newboxsize
    # check if atoms go over newboxsize
    lim_match = {min: 'lo', max: 'hi'}
    merged_limits = limitsTopology(merged_topology)
    for limit in merged_limits:
        for coord in merged_limits[limit]:
            under_test = merged_limits[limit][coord]
            if (limit(under_test,
                      newboxsize[coord+lim_match[limit]]) == under_test) and (
                      # when both are the same
                      under_test != newboxsize[coord+lim_match[limit]]):
                warn("Atoms have coordinates outside new box. \
                       Increase boxsize.", UserWarning)
    # cleanup
    del merged_topology['topologycounts'][None]
    return merged_topology


def writeTopology(topology, filename):
    """
    Writes topology to a file in the LAMMPS data format.

    Parameters
    ----------
    topology : nested dicts
        Dictionary in the style of atoms.getAllAtomData().

    filename : str
        Name of the datafile to be write topology to.
    """

    names = topology['atomnames']
    badis = {'bonddata': "Bonds",
             'angledata': "Angles",
             'dihedraldata': "Dihedrals",
             'improperdata': "Impropers"}

    with open(filename, "w") as f:
        # LAMMPS ignores the first line. Let's print some useful data. Max 254
        # characters.
        date = (datetime.now()
                        .astimezone()
                        .replace(microsecond=0)
                        .isoformat(' ')
                )
        firstline = (f"LAMMPS data file merged using {_name}.{_this_file} "
                     f"version {__version__} "
                     f"on {date}"
                     )[:254]

        f.write(firstline)
        # .write() does not do \n automatically:
        f.write("\n")
        # Header
        for propert in topology['topologycounts']:
            f.write(f"{topology['topologycounts'][propert]} {propert}\n")
        # Box size
        box = topology['boxsize']
        box_order = ['xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi']
        for low, high in zip(box_order[::2], box_order[1::2]):
            f.write(f"{box[low]} {box[high]} {low} {high}\n")

        # Pair/Bond/Angle/Dihedral/Improper Coeffs
        for coeff in coeffs_equiv:
            # Write only if it's > 0
            if topology['topologycounts'][atoms.g_coeffs[coeff]] > 0:
                # start section, separate from previous
                f.write("\n")
                # using g_coeffs because lammps uses pair types and atom types
                # using always the same designation would make things too easy.
                f.write(coeff)
                f.write("\n")
                # Empty line
                f.write("\n")
                section = topology[coeffs_equiv[coeff]]
                for line in section:
                    # id
                    params = [str(a) for a in section[line][:-1]]
                    f.write(f"{str(line)}")
                    f.write(" ")
                    f.write(f"{' '.join(params)}")
                    f.write(f" # {section[line][-1].lstrip()}")
                    f.write("\n")

        # Masses
        if len(topology['masses']) > 0:
            # start section, separate from previous
            f.write("\n")
            f.write("Masses")
            f.write("\n")
            # Empty line
            f.write("\n")
            for atomtype in topology['atomnames']:
                # id
                f.write(f"{atomtype}")
                f.write(" ")
                # mass
                f.write(f"{topology['masses'][names[atomtype]]:.4f}")
                # atom type str
                f.write(f" # {topology['atomnames'][atomtype]}")
                f.write("\n")

        # Atom data
        if len(topology['atomdata']) > 0:
            # start section, separate from previous
            f.write("\n")
            f.write("Atoms")
            f.write("\n")
            # Empty line
            f.write("\n")
            for atom in topology['atomdata']:
                f.write(f"{atom}")
                f.write(" ")
                f.write(f"{topology['atomdata'][atom]['mol']}")
                f.write(" ")
                f.write(f"{topology['atomdata'][atom]['type']}")
                f.write(" ")
                # 4 decimal places for charge
                f.write(f"{topology['atomdata'][atom]['charge']:.4f}")
                f.write(" ")
                # 6 decimal places for coordinates
                f.write(f"{topology['atomdata'][atom]['x']:.6f}")
                f.write(" ")
                f.write(f"{topology['atomdata'][atom]['y']:.6f}")
                f.write(" ")
                f.write(f"{topology['atomdata'][atom]['z']:.6f}")
                # atom type str
                f.write(f" # {names[topology['atomdata'][atom]['type']]}")
                f.write("\n")

        # Bond/Angle/Dihedral/Improper data
        f.write(_getStringBADI(topology, badis))

    return


###############################################################################
#                             Protected Functions                             #
###############################################################################


def _getStringBADI(topology, badis):
    # writeTopology helper function, takes topology dict and badis dict,
    # returns string of all the lines in all BADIS ready to be written to file.
    string = ''
    for badi in badis:
        if len(topology[badi]) > 0:
            string += str('\n' + badis[badi]) + '\n\n'
            for line in topology[badi]:
                string += str(line) + ' ' + ' '.join(
                    [str(a) for a in topology[badi][line]]) + '\n'
    return string


def _shiftKey(topology_property: dict, shiftkey: int):
    # adds shiftkey to dictionarty key (int)
    return {k + shiftkey: v for k, v in topology_property.items()}


def _shiftDictOfLists(topology_property,
                      shiftkey,
                      shiftlistid,
                      shiftlistelements,
                      molid):
    """
    shifts dict of lists.
    {key + shiftkey:
        [listid + shiftlistid,
        el + shiftlistelements,
        el + shiftlistelements, etc]}
    """

    if shiftlistid is not None:
        new_topology_property = {}
        for item in topology_property:
            new_list = topology_property[item].copy()
            new_list[0] += shiftlistid
            new_list[1:] = (x + shiftlistelements for x in new_list[1:])
            new_topology_property[item + shiftkey] = new_list
        return new_topology_property
        #  return {k + shiftkey: [a for a in v] \
        #  for k, v in topology_property.items()}
    else:
        return _shiftKey(topology_property, shiftkey)


def _shiftDictOfDicts(topology_property,
                      shiftkey,
                      targetinnerkey,
                      shiftinnervalue,
                      molid):
    """
    shifts dict of lists.
    {key + shiftkey:
        {unrelated_key: unrelated_value,
        targetinnerkey: original_value + shiftinnervalue,
        unrelated_key: unrelated_value}}
    """

    new_topology_property = _shiftKey(topology_property, shiftkey)
    new_topology_property = _setMolid(new_topology_property, molid)
    new_topology_property = _addToAtomProperty(new_topology_property,
                                               targetinnerkey,
                                               shiftinnervalue)
    return new_topology_property


def _addToAtomProperty(topology, propert, value):
    # adds value to every topology[:][propert]
    for atom in topology:
        topology[atom][propert] += value
    return topology


def _setMolid(topology, molid):
    # changes molid
    for atom in topology:
        topology[atom]['mol'] = molid
    return topology


def _combineDicts(a, b, op=operator.add):
    # dark magic. From: https://stackoverflow.com/a/11012181
    return dict(list(a.items()) + list(b.items()) +
                [(k, op(a[k], b[k])) for k in set(b) & set(a)])


def _mergeDicts(a, b):
    # for mass, warns about similar keys with different values.
    error = 'Trying to merge two dictionaries that share a key with'\
            ' different values.'
    result = deepcopy(a)
    try:
        for key in b:
            if key in a.keys() and a[key] != b[key]:
                raise ValueExists
            else:
                result[key] = b[key]
    except ValueExists as e:
        message = f"{error} Offending key: {key}"
        print(message)
        raise ValueExists(message) from e

    return result


def _strToFunction(settings):
    # changes min/max (str) to min/max (function)
    for key in settings:
        func = settings[key]['minormax']
        if func == 'min':
            settings[key]['minormax'] = min
        elif func == 'max':
            settings[key]['minormax'] = max
        else:
            func = 'None' if func == '' else func
            raise ValueError(
                f"The key 'minormax' accepts only 'min' or 'max', not {func}")
    return settings


def _openYaml(ifile):
    # to read settings file
    with open(ifile, 'r') as f:
        settings = full_load(f)
    return _strToFunction(settings['inputs']), settings['output']


def _limit(function, coordinate, topology):
    # gets min/max coordinates in topology
    return topology['atomdata'][function(
        topology['atomdata'], key=lambda atom:
            topology['atomdata'][atom][coordinate])][coordinate]


def main(args=argv):
    print("")
    inputs, outputsettings = readInputFile(args)
    # read files
    topologies = readTopologies(inputs)
    # find minimums and maximums
    limits = limitsAllTopologies(topologies)
    # shift mins/maxs to defined values
    new_topologies = shiftTopologies(topologies, limits, inputs)
    #  create new merged topology
    merged_topology = mergeTopologies(
        new_topologies, outputsettings['boxsize'])
    # write merged topology to datafile
    writeTopology(merged_topology, outputsettings['filename'])
    print("File created Successfully")
    print("")
