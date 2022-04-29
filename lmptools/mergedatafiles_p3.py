#!/usr/bin/python3
#  import lmptools.readfiles as rdfl
import lmptools.atoms as atoms
#  import lmptools.commondata_p3 as cdp3
from yaml import full_load
from sys import argv


class MissingSettingsFile(IOError):
    pass


class TooManyArguments(IOError):
    pass


args = argv


def readInputFile(args):
    try:
        if len(args) < 2:
            for filename in ['merge.yaml', 'merge.yml']:
                try:
                    return _openYaml(filename)

                except IOError:
                    continue
            else:
                raise MissingSettingsFile
        elif len(args) > 2:
            raise TooManyArguments
        else:
            filename = args[1]
            return _openYaml(filename)

    except MissingSettingsFile:
        # TODO expand instructions
        _myExit("Missing settings file: merge.yaml or merge.yml", 3)

    except TooManyArguments:
        # TODO expand instructions
        _myExit("This script takes only one argument, the settings file.", 4)


def readTopologies(files):
    topologies = {}
    for file in files:
        topologies[file] = atoms.getAllAtomData(file['filename'])
        # check max Z
        #  a[max(a, key=lambda v: a[v]['z'])]['z']
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
        topology, key=lambda atom:
            topology['atomdata'][atom][coordinate])][coordinate]


if __name__ == '__main__':
    files = readInputFile()  # pragma: no cover
    # read files
    topologies = readTopologies(files)
    # write new datafile

    # TODO: consider pytest-console-scripts tests
