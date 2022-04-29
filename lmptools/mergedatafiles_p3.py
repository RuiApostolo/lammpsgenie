#!/usr/bin/python3
#  import lmptools.readfiles as rdfl
#  import lmptools.atoms as atoms
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
                    settings = _openYaml(filename)
                    break
                except IOError:
                    continue
            else:
                raise MissingSettingsFile
        elif len(args) > 2:
            raise TooManyArguments
        else:
            filename = args[1]
            settings = _openYaml(filename)

    except MissingSettingsFile:
        # TODO expand instructions
        _myExit("Missing settings file: merge.yaml or merge.yml", 3)

    except TooManyArguments:
        _myExit("This script takes only one argument, the settings file.", 4)

    return settings


def _myExit(message, code):
    print(message)
    exit(code)


def _openYaml(ifile):
    print(ifile)
    with open(ifile, 'r') as f:
        settings = full_load(f)
    print(settings)
    return settings


if __name__ == '__main__':
    readInputFile()  # pragma: no cover
    # TODO: consider pytest-console-scripts tests
