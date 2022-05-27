#!/usr/bin/env python3
from pkgutil import get_data

iron_50 = ['Fe2O3_50_down.lammps', 'Fe2O3_50_up.lammps']
iron_100 = ['Fe2O3_100_down.lammps', 'Fe2O3_100_up.lammps']


def save_iron_50():
    return _save_files(iron_50)


def save_iron_100():
    return _save_files(iron_100)


def save_iron_all():
    save_iron_50()
    save_iron_100()


def _save_files(flist):
    """
    Helper function that saves files.
    """

    for file in flist:
        data = get_data(__name__, "data/"+file).decode('utf8')
        with open(file, "w") as f:
            f.write(data)
