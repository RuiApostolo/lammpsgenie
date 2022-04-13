import pytest
from conftest import zipRefs, ref_data_files
import lmptools.atoms as atoms
import lmptools.commondata_p3 as cdp3


def test_getNatoms_fromDummy(dummy_dump):
    assert cdp3.getNatoms(dummy_dump) == 16889


@pytest.mark.parametrize("natoms", [2400])
def test_getNatoms_fromFile(dumpfilelines, natoms):
    assert cdp3.getNatoms(dumpfilelines) == natoms


class TestAtomData:
    ref_atom_types = [{1: 'SCP', 2: 'SCS'},
                      {1: 'CHK', 2: 'COK', 3: 'HCK', 4: 'OCK'},
                      {1: 'COO', 2: 'HCO', 3: 'OCO'},
                      {1: 'CHE', 2: 'COE', 3: 'HCE', 4: 'HOE', 5: 'OHE'}]

    ref_atom_data = [{'mol': 1,
                      'type': 'SCP',
                      'charge': 0.0,
                      'x': 36.114258,
                      'y': 28.328382,
                      'z': 34.113575},
                     {'mol': 1,
                      'type': 'CHK',
                      'charge': -0.25,
                      'x': -2.1294,
                      'y': 0.4027,
                      'z': 0.0},
                     {'mol': 1,
                      'type': 'COO',
                      'charge': 0.135,
                      'x': 1.00761,
                      'y': -0.2128,
                      'z': -0.182370},
                     {'mol': 1,
                      'type': 'CHE',
                      'charge': -0.21,
                      'x': 1.092820,
                      'y': -0.05313,
                      'z': 0.09623},
                     ]

    ref_masses = [{'SCP': 15.035, 'SCS': 14.027},
                  {'CHK': 12.011, 'COK': 12.011, 'HCK': 1.008, 'OCK': 15.999},
                  {'COO': 12.011, 'HCO': 1.008, 'OCO': 15.999},
                  {'CHE': 12.011, 'COE': 12.011, 'HCE': 1.008,
                   'HOE': 1.008, 'OHE': 15.999},
                  ]

    ref_boxsizes = [{'xlo': 0.0, 'xhi': 45.0,
                     'ylo': 0.0, 'yhi': 45.0,
                     'zlo': 0.0, 'zhi': 45.0},
                    {'xlo': -3.4836, 'xhi': 1.5164,
                     'ylo': -2.0825, 'yhi': 2.9175,
                     'zlo': -2.5, 'zhi': 2.5},
                    {'xlo': -0.74882, 'xhi': 4.25118,
                     'ylo': -2.405085, 'yhi': 2.594915,
                     'zlo': -2.24183, 'zhi': 2.75817},
                    {'xlo': -0.416785, 'xhi': 4.583215,
                     'ylo': -2.552375, 'yhi': 2.447625,
                     'zlo': -2.32406, 'zhi': 2.67594},
                    ]

    @pytest.mark.parametrize("datafile, types",
                             zipRefs(ref_data_files(),
                                     ref_atom_types))
    def test_getAtomType(self, datafile, types):
        assert atoms.getAtomType(datafile) == types

    @pytest.mark.parametrize("datafile, atomdata, masses",
                             zipRefs(ref_data_files(),
                                     ref_atom_data,
                                     ref_masses))
    def test_getAtomData(self, datafile, atomdata, masses):
        assert atoms.getAtomData(datafile)[0][1] == atomdata
        assert atoms.getAtomData(datafile)[1] == masses

    @pytest.mark.parametrize("datafile, atomdata, masses, boxsizes",
                             zipRefs(ref_data_files(),
                                     ref_atom_data,
                                     ref_masses,
                                     ref_boxsizes))
    def test_getAllAtomData(self, datafile, atomdata, masses, boxsizes):
        assert atoms.getAllAtomData(datafile)[0][1] == atomdata
        assert atoms.getAllAtomData(datafile)[5] == masses
        assert atoms.getAllAtomData(datafile)[6] == boxsizes
