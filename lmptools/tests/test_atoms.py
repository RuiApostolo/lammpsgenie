import pytest
import lmptools.atoms as atoms
import lmptools.commondata_p3 as cdp3


def test_getNatoms_fromDummy(dummy_dump):
    assert cdp3.getNatoms(dummy_dump) == 16889


@pytest.mark.parametrize("natoms", [2400])
def test_getNatoms_fromFile(dumpfilelines, natoms):
    assert cdp3.getNatoms(dumpfilelines) == natoms


@pytest.mark.parametrize(
    "datafile, atomtypes", [
         ("tests/uadodecane.data", {1: 'SCP', 2: 'SCS'})
    ])
def test_getAtomType(datafile, atomtypes):
    assert atoms.getAtomType("tests/uadodecane.data") == atomtypes


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
    assert atoms.getAtomData(filename)[0][1] == atomdata
    assert atoms.getAtomData(filename)[1] == masses


@pytest.mark.parametrize(
    "filename, atomdata, masses, boxsizes", [
     ("tests/uadodecane.data",
      {'mol': 1,
       'type': 'SCP',
       'charge': 0.0,
       'x': 36.114258,
       'y': 28.328382,
       'z': 34.113575},
      {'SCP': 15.035,
       'SCS': 14.027},
      {'xlo': 0.0,
       'xhi': 45.0,
       'ylo': 0.0,
       'yhi': 45.0,
       'zlo': 0.0,
       'zhi': 45.0}
      )
    ])
def test_getAllAtomData(filename, atomdata, masses, boxsizes):
    assert atoms.getAllAtomData(filename)[0][1] == atomdata
    assert atoms.getAllAtomData(filename)[5] == masses
    assert atoms.getAllAtomData(filename)[6] == boxsizes
