#!/usr/bin/python3
#  import lmptools.readfiles as rdfl
import lmptools.atoms as atoms
#  import lmptools.commondata_p3 as cdp3
from yaml import full_load
from sys import argv, exit


class MissingSettingsFile(IOError):
    pass


class TooManyArguments(IOError):
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

    except MissingSettingsFile:
        # TODO expand instructions
        _myExit("Missing settings file: merge.yaml or merge.yml", 3)

    except TooManyArguments:
        # TODO expand instructions
        _myExit("This script takes only one argument, the settings file.", 4)


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

    # what to add to current coordinates
    delta = settings['value'] - \
        settings['minormax'](limits[min][axis], limits[max][axis])
    for atom in topology['atomdata']:
        topology['atomdata'][atom][axis] += delta
    return topology


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
    return _strToFunction(settings)


def _limit(function, coordinate, topology):
    return topology['atomdata'][function(
        topology['atomdata'], key=lambda atom:
            topology['atomdata'][atom][coordinate])][coordinate]


if __name__ == '__main__':  # pragma: no cover
    settings = readInputFile()
    # read files
    topologies = readTopologies(settings)
    # find minimums and maximums
    limits = limitsAllTopologies(topologies)
    # shift mins/maxs to defined values
    new_topologies = shiftTopologies(topologies, limits, settings)
    # write new datafile

    # TODO: consider pytest-console-scripts tests
