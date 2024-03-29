import pytest
import lammpsgenie.commondata as cdt
import lammpsgenie.atoms as atoms


def test_getDumpTSRange_fromDummy(dummy_dump):
    assert cdt.getDumpTSRange(dummy_dump) == [1000000, 2000000]


@pytest.mark.parametrize(
    "timesteps", [
     [10000, 20000, 30000, 40000, 50000, 60000]
     ])
def test_getDumpTSRange_fromFile(dumpfilelines, timesteps):
    assert cdt.getDumpTSRange(dumpfilelines) == timesteps


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
    traj = cdt.readTS(dumpfilelines,
                      framenumber,
                      atoms.getAtomType("tests/uadodecane.data"),
                      9
                      )
    frame = cdt.getDumpTSRange(dumpfilelines)[framenumber-1]
    assert len(traj.keys()) == 1
    assert len(traj[frame]['atom'].keys()) == 2400
    assert traj[frame]['boxsize'] == boxsize
    assert traj[frame]['atom'][atomnumber] == properties


@pytest.mark.parametrize(
    "start, stop, numberframes", [
        ('first', 'last', 6),
        ('first', 1, 1),
        ('first', 3, 3),
        (4, 'last', 3),
        (6, 'last', 1),
        (2, 2, 1),
        (3, 2, 0),
        (10, 2, 0),
    ])
def test_getTrajTSRange(dumpfilelines,
                        start,
                        stop,
                        numberframes):
    traj = {}
    for frame in range(1, 6 + 1):
        traj[frame] = cdt.readTS(
                          dumpfilelines,
                          frame,
                          atoms.getAtomType("tests/uadodecane.data"),
                          9
                          )
    slice_traj = cdt.getTrajTSRange(traj, start, stop)
    assert len(slice_traj) == numberframes
