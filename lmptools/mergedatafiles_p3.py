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
    settings : List of dicts
        Each dict should contain:
        'filename': str
        'minormax': str
        'value': floa

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

    return atoms.getAllAtomData(file['filename'])


def readTopologies(files):
    """
    Reads Topology from settings dict.

    Parameters
    ----------
    files : list of dicts
        Settings, from readInputFile().

    Returns
    -------
    topology : nested dicts
        dictionary of the form:
        'filename': Dictionary in the style of atoms.getAllAtomData()
    """

    topologies = {}
    for file in files:
        topologies[file['filename']] = readTopology(file)
    return topologies


def _myExit(message, code):
    print(message)
    exit(code)


def _openYaml(ifile):
    with open(ifile, 'r') as f:
        settings = full_load(f)
    return settings


def _limit(function, coordinate, topology):
    return topology['atomdata'][function(
        topology['atomdata'], key=lambda atom:
            topology['atomdata'][atom][coordinate])][coordinate]


if __name__ == '__main__':  # pragma: no cover
    files = readInputFile()
    # read files
    topologies = readTopologies(files)
    # write new datafile

    # TODO: consider pytest-console-scripts tests
