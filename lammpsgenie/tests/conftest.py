import pytest
import lammpsgenie.readfiles as rdfl


@pytest.fixture
def emptylist():
    return []


@pytest.fixture
def mock_path(monkeypatch):
    monkeypatch.chdir("tests")


@pytest.fixture
def dummy_dump():
    return [
         "ITEM: TIMESTEP",
         "1000000",
         "ITEM: NUMBER OF ATOMS",
         "16889",
         "ITEM: BOX BOUNDS pp pp ff",
         "0.0000000000000000e+00 5.5088000000000001e+01",
         "0.0000000000000000e+00 5.0380000000000003e+01",
         "-2.0000000000000000e+01 1.9000000000000000e+02",
         "ITEM: TIMESTEP",
         "2000000",
         "ITEM: NUMBER OF ATOMS",
         "16889",
         "ITEM: BOX BOUNDS pp pp ff",
         "0.0000000000000000e+00 5.5088000000000001e+01",
         "0.0000000000000000e+00 5.0380000000000003e+01",
         "-2.0000000000000000e+01 1.9000000000000000e+02"
         ]


@pytest.fixture
def dumpfilelines():
    tdump = rdfl.readAll("tests/dump.uadodecane.lammpstrj")
    return tdump


def ref_data_files():
    datafiles = (
        "ketene.lammps",
        "oxirene.lammps",
        "ethynol.lammps",
        "oxirene_bare.lammps"
                 )
    return _add_testd(datafiles)


def ref_data_files_large():
    return _add_testd(["uadodecane.data"]) + ref_data_files()


def ref_iron_files():
    ironfiles = (
        "data/Fe2O3_50_down.lammps",
        "data/Fe2O3_50_up.lammps",
        "data/Fe2O3_100_down.lammps",
        "data/Fe2O3_100_up.lammps",
    )
    return [a for a in ironfiles]


def ref_all_data_fs():
    return ref_data_files_large() + ref_iron_files()


def _add_testd(files: list):
    return ["tests/" + a for a in files]


def zipRefs(*args):
    return list(zip(*args))
