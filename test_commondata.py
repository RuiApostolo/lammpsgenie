import pytest
import commondata_p3 as cdp3


@pytest.mark.parametrize(
    "filename, lines", [
        ("Fe2O3_50_down.lammps", 14120),
        ("Fe2O3_50_up.lammps", 14120),
        ("Fe2O3_100_down.lammps", 56360),
        ("Fe2O3_100_up.lammps", 56360)
    ])
def test_readAll_fromParameters(filename, lines):
    assert len(cdp3.readAll(filename)) == lines


@pytest.mark.parametrize(
    "filename, nlines", [
        ("dump.uadodecane.lammpstrj", 14454)
    ])
def test_readAll_fromFile(filename, nlines):
    assert len(cdp3.readAll(filename)) == nlines


@pytest.mark.parametrize(
    "filename, lines", [
        ("Fe2O3_50_down.lammps.gz", 14120),
        ("Fe2O3_50_up.lammps.gz", 14120),
        ("Fe2O3_100_down.lammps.gz", 56360),
        ("Fe2O3_100_up.lammps.gz", 56360)
    ])
def test_readAllGzip(filename, lines):
    assert len(cdp3.readAllGzip(filename)) == lines


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


def test_getNatoms_fromDummy(dummy_dump):
    assert cdp3.getNatoms(dummy_dump) == 16889


@pytest.fixture
def dumpfilelines():
    tdump = cdp3.readAll("dump.uadodecane.lammpstrj")
    return tdump


@pytest.mark.parametrize("natoms", [2400])
def test_getNatoms_fromFile(dumpfilelines, natoms):
    assert cdp3.getNatoms(dumpfilelines) == natoms


def test_getTSrange_fromDummy(dummy_dump):
    assert cdp3.getTSrange(dummy_dump) == [1000000, 2000000]


@pytest.mark.parametrize(
    "timesteps", [
     [10000, 20000, 30000, 40000, 50000, 60000]
     ])
def test_getTSrange_fromFile(dumpfilelines, timesteps):
    assert cdp3.getTSrange(dumpfilelines) == timesteps


@pytest.mark.parametrize(
    "datafile, atomtypes", [
         ("uadodecane.data", {'1': 'SCP', '2': 'SCS'})
    ])
def test_getAtomType(datafile, atomtypes):
    assert cdp3.getAtomType("uadodecane.data") == atomtypes


#  def test_readTS_fromFile(dumpfilelines):
#      traj = cdp3.readTS(dumpfilelines,
#                         9,
#                         cdp3.getNatoms(dumpfilelines),
#                         len(cdp3.getTSrange(dumpfilelines)),
#                         )
#      assert len(traj.keys()) == 6

# TODO: getTS
