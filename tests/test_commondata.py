import pytest
import src.commondata_p3 as cdp3
#  import lmpdttrj.commondata_p3 as cdp3


@pytest.mark.parametrize(
    "filename, lines", [
        ("tests/Fe2O3_50_down.lammps", 14120),
        ("tests/Fe2O3_50_up.lammps", 14120),
        ("tests/Fe2O3_100_down.lammps", 56360),
        ("tests/Fe2O3_100_up.lammps", 56360)
    ])
def test_readAll_fromParameters(filename, lines):
    assert len(cdp3.readAll(filename)) == lines


@pytest.mark.parametrize(
    "filename, nlines", [
        ("tests/dump.uadodecane.lammpstrj", 14454)
    ])
def test_readAll_fromFile(filename, nlines):
    assert len(cdp3.readAll(filename)) == nlines


@pytest.mark.parametrize(
    "filename, lines", [
        ("tests/Fe2O3_50_down.lammps.gz", 14120),
        ("tests/Fe2O3_50_up.lammps.gz", 14120),
        ("tests/Fe2O3_100_down.lammps.gz", 56360),
        ("tests/Fe2O3_100_up.lammps.gz", 56360)
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
    tdump = cdp3.readAll("tests/dump.uadodecane.lammpstrj")
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
         ("tests/uadodecane.data", {1: 'SCP', 2: 'SCS'})
    ])
def test_getAtomType(datafile, atomtypes):
    assert cdp3.getAtomType("tests/uadodecane.data") == atomtypes


@pytest.mark.parametrize(
    "framenumber, atomnumber, boxsize, properties", [
        (1,
         1,
         (4.5000000000000000e+01 - 0.0000000000000000e+00,
          4.5000000000000000e+01 - 0.0000000000000000e+00,
          4.5000000000000000e+01 - 0.0000000000000000e+00),
         {'type': 'SCP',
          'x': 36.5534,
          'y': 30.3231,
          'z': 35.2407,
          'vx': -0.00334412,
          'vy': -0.0064402,
          'vz': -0.000754788}
         ),
        (2,
         13,
         (4.3547785489204998e+01 - 1.4522145107952262e+00,
          4.3547785489204998e+01 - 1.4522145107952262e+00,
          4.3547785489204998e+01 - 1.4522145107952262e+00),
         {'type': 'SCP',
          'x': 5.8559,
          'y': 38.4064,
          'z': 41.7078,
          'vx': -0.00277012,
          'vy': 0.00170807,
          'vz': -0.000710621}
         ),
        (3,
         117,
         (4.3382757517072228e+01 - 1.6172424829280736e+00,
          4.3382757517072228e+01 - 1.6172424829280736e+00,
          4.3382757517072228e+01 - 1.6172424829280736e+00),
         {'type': 'SCS',
          'x': 41.6643,
          'y': 34.7271,
          'z': 29.4161,
          'vx': 0.00171461,
          'vy': -0.00668507,
          'vz': 0.0018407}
         ),
        (4,
         1311,
         (4.3436667898165098e+01 - 1.5633321018358153e+00,
          4.3436667898165098e+01 - 1.5633321018358153e+00,
          4.3436667898165098e+01 - 1.5633321018358153e+00),
         {'type': 'SCS',
          'x': 16.664,
          'y': 36.9752,
          'z': 32.6149,
          'vx': 0.00143178,
          'vy': -0.00434612,
          'vz': 0.00847797}
         ),
        (5,
         2155,
         (4.3335642577939524e+01 - 1.6643574220614781e+00,
          4.3335642577939524e+01 - 1.6643574220614781e+00,
          4.3335642577939524e+01 - 1.6643574220614781e+00),
         {'type': 'SCS',
          'x': 39.211,
          'y': 15.9501,
          'z': 4.87282,
          'vx': -0.00469417,
          'vy': -0.00120414,
          'vz': 0.00573282}
         ),
        (6,
         2400,
         (4.3358792843046842e+01 - 1.6412071569538185e+00,
          4.3358792843046842e+01 - 1.6412071569538185e+00,
          4.3358792843046842e+01 - 1.6412071569538185e+00),
         {'type': 'SCP',
          'x': 38.8128,
          'y': 6.34872,
          'z': 41.3552,
          'vx': 0.006193,
          'vy': -0.000873072,
          'vz': -0.0010738}
         )
    ])
def test_readTS_fromFile(dumpfilelines,
                         framenumber,
                         atomnumber,
                         boxsize,
                         properties):
    traj = cdp3.readTS(dumpfilelines,
                       framenumber,
                       cdp3.getAtomType("tests/uadodecane.data"),
                       9
                       )
    frame = cdp3.getTSrange(dumpfilelines)[framenumber-1]
    assert len(traj.keys()) == 1
    assert len(traj[frame]['atom'].keys()) == 2400
    assert traj[frame]['boxsize'] == boxsize
    assert traj[frame]['atom'][atomnumber] == properties

# TODO: getTS


@pytest.mark.parametrize(
    "filename, atomdata, masses", [
     ("tests/uadodecane.data",
      {'mol': 1,
       'type': 'SCP',
       'charge': 0.0,
       'x': 36.114258,
       'y': 28.328382,
       'z': 34.113575},
      {'SCP': 15.035,
       'SCS': 14.027}
      )
    ])
def test_getAtomData(filename, atomdata, masses):
    assert cdp3.getAtomData(filename)[0][1] == atomdata
    assert cdp3.getAtomData(filename)[1] == masses
