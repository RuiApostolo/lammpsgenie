#!/usr/bin/python3
#  import lmptools.readfiles as rdfl
import lmptools.atoms as atoms
#  import lmptools.commondata_p3 as cdp3
from yaml import full_load
from sys import argv, exit
from copy import deepcopy
import operator


class MissingSettingsFile(IOError):
    pass


class TooManyArguments(IOError):
    pass


class ValueExists(ValueError):
    pass


args = argv


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
        # TODO expand instructions
        message = "Missing settings file: merge.yaml or merge.yml"
        raise MissingSettingsFile(message) from exception
        _myExit(message, 3)

    except TooManyArguments as exception:
        # TODO expand instructions
        message = "This script takes only one argument, the settings file."
        raise MissingSettingsFile(message) from exception
        _myExit(message, 4)


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
    Reads Topology from settings dict.

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
        With shifted coordinates
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


def mergeTopologies(topologies):
    new_topology = {}
    property_pairs = {
        'pairtypes': 'atom types',
        'bondtypes': 'bond types',
        'angletypes': 'angle types',
        'dihedraltypes': 'dihedral types',
        'impropertypes': 'improper types',
        'atomdata': 'atoms',
        'bonddata': 'bonds',
        'angledata': 'angles',
        'dihedraldata': 'dihedrals',
        'improperdata': 'impropers',
    }
    for topology in topologies:
        if len(new_topology) == 0:
            new_topology = deepcopy(topologies[topology])
        else:
            shift_values = new_topology['topologycounts']
            for propert in property_pairs:
                new_topology[propert] = _shiftBy(
                    topology[propert],
                    shift_values[property_pairs[propert]])

            # add topology counts
            new_topology['topologycounts'] = _combineDicts(
                shift_values, topologies[topology]['topologycounts'])
    # limits/boxsize
    return new_topology


def _shiftBy(topology_property, count):
    return {k + count: v for k, v in topology_property.items()}


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
        raise ValueExists(message) from e
        _myExit(message, 5)

    return result


def _myExit(message, code):
    print(message)
    exit(code)


def _strToFunction(settings):
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
    with open(ifile, 'r') as f:
        settings = full_load(f)
    return _strToFunction(settings['inputs']), settings['output']


def _limit(function, coordinate, topology):
    return topology['atomdata'][function(
        topology['atomdata'], key=lambda atom:
            topology['atomdata'][atom][coordinate])][coordinate]


if __name__ == '__main__':  # pragma: no cover
    settings, outputfile = readInputFile()
    # read files
    topologies = readTopologies(settings)
    # find minimums and maximums
    limits = limitsAllTopologies(topologies)
    # shift mins/maxs to defined values
    new_topologies = shiftTopologies(topologies, limits, settings)
    #  create new merged topology
    merged_topology = mergeTopologies(new_topologies)
    # write merged topology to datafile
    #  writeTopology(outputfile)
    # TODO: consider pytest-console-scripts tests
