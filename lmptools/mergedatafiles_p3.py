#!/usr/bin/python3
#  import lmptools.readfiles as rdfl
#  import lmptools.atoms as atoms
#  import lmptools.commondata_p3 as cdp3
import yaml
from sys import argv


class MissingSettingsFile(IOError):
    pass


class TooManyArguments(IOError):
    pass


try:
    if len(argv) < 2:
        for ifile in ['merge.yaml', 'merge.yml']:
            try:
                with open(ifile, 'r') as f:
                    settings = yaml.full_load(f)
            except IOError:
                continue
        else:
            raise MissingSettingsFile
    elif len(argv) > 2:
        raise TooManyArguments
    else:
        with open(ifile, 'r') as f:
            settings = yaml.full_load(f)

except MissingSettingsFile:
    # TODO expand instructions
    print('Missing settings file: merge.yaml or merge.yml')
    exit(3)

except TooManyArguments:
    print("This script takes only one argument, the settings file.")
    exit(4)
