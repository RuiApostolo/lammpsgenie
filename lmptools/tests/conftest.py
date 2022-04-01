import pytest
import lmptools.readfiles as rdfl


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
