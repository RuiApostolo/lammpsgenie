import pytest
from conftest import zipRefs, ref_data_files, ref_data_files_large
import lmptools.atoms as atoms
import lmptools.commondata_p3 as cdp3


def test_getNatoms_fromDummy(dummy_dump):
    assert cdp3.getNatoms(dummy_dump) == 16889


@pytest.mark.parametrize("natoms", [2400])
def test_getNatoms_fromFile(dumpfilelines, natoms):
    assert cdp3.getNatoms(dumpfilelines) == natoms


class TestAtomDataLarge:
    ref_atom_types = [
        # ua_dodecane
        {1: 'SCP', 2: 'SCS'},
        # ketene
        {1: 'CHK', 2: 'COK', 3: 'HCK', 4: 'OCK'},
        # oxirene
        {1: 'COO', 2: 'HCO', 3: 'OCO'},
        # ethynol
        {1: 'CHE', 2: 'COE', 3: 'HCE', 4: 'HOE', 5: 'OHE'},
        # oxirene_bare
        {1: 'COO', 2: 'HCO', 3: 'OCO'},
                      ]

    ref_atom_data = [
        # ua_dodecane
        {'mol': 1,
         'type': 'SCP',
         'charge': 0.0,
         'x': 36.114258,
         'y': 28.328382,
         'z': 34.113575},
        # ketene
        {'mol': 1,
         'type': 'CHK',
         'charge': -0.25,
         'x': -2.1294,
         'y': 0.4027,
         'z': 0.0},
        # oxirene
        {'mol': 1,
         'type': 'COO',
         'charge': 0.135,
         'x': 1.00761,
         'y': -0.2128,
         'z': -0.182370},
        # ethynol
        {'mol': 1,
         'type': 'CHE',
         'charge': -0.21,
         'x': 1.092820,
         'y': -0.05313,
         'z': 0.09623},
        # oxirene_bare
        {'mol': 1,
         'type': 'COO',
         'charge': 0.135,
         'x': 1.00761,
         'y': -0.2128,
         'z': -0.182370},
                     ]

    ref_masses = [
        # ua_dodecane
        {'SCP': 15.035, 'SCS': 14.027},
        # ketene
        {'CHK': 12.011, 'COK': 12.011, 'HCK': 1.008, 'OCK': 15.999},
        # oxirene
        {'COO': 12.011, 'HCO': 1.008, 'OCO': 15.999},
        # ethynol
        {'CHE': 12.011, 'COE': 12.011, 'HCE': 1.008,
         'HOE': 1.008, 'OHE': 15.999},
        # oxirene_bare
        {'COO': 12.011, 'HCO': 1.008, 'OCO': 15.999},
                  ]

    ref_boxsizes = [
        # ua_dodecane
        {'xlo': 0.0, 'xhi': 45.0,
         'ylo': 0.0, 'yhi': 45.0,
         'zlo': 0.0, 'zhi': 45.0},
        # ketene
        {'xlo': -3.4836, 'xhi': 1.5164,
         'ylo': -2.0825, 'yhi': 2.9175,
         'zlo': -2.5, 'zhi': 2.5},
        # oxirene
        {'xlo': -0.74882, 'xhi': 4.25118,
         'ylo': -2.405085, 'yhi': 2.594915,
         'zlo': -2.24183, 'zhi': 2.75817},
        # ethynol
        {'xlo': -0.416785, 'xhi': 4.583215,
         'ylo': -2.552375, 'yhi': 2.447625,
         'zlo': -2.32406, 'zhi': 2.67594},
        # oxirene_bare
        {'xlo': -0.74882, 'xhi': 4.25118,
         'ylo': -2.405085, 'yhi': 2.594915,
         'zlo': -2.24183, 'zhi': 2.75817},
                    ]

    ref_pair_coeffs = [
        # ua_dodecane
        {1: [0.175, 3.905, ' 1  SCP 10 '],
         2: [0.118, 3.905, ' 2  SCS 13 ']},
        # ketene
        {1: [0.086, 3.3, ' 1  CHK 900/47'],
         2: [0.086, 3.3, ' 2  COK 904/110'],
         3: [0.03, 2.42, ' 3  HCK 899/46'],
         4: [0.21, 2.96, ' 4  OCK 905/4']},
        # oxirene
        {1: [0.076, 3.75, ' 1  COO 87/47'],
         2: [0.03, 2.42, ' 2  HCO 89/46'],
         3: [0.14, 2.9, ' 3  OCO 122/20']},
        # ethynol
        {1: [0.086, 3.3, ' 1  CHE 755/19'],
         2: [0.1, 3.3, ' 2  COE 759/19'],
         3: [0.015, 2.42, ' 3  HCE 756/46'],
         4: [0.0, 0.0, ' 4  HOE 97/7'],
         5: [0.17, 3.12, ' 5  OHE 96/5']},
        # oxirene_bare
        {1: [' 1  COO 87/47'],
         2: [' 2  HCO 89/46'],
         3: [' 3  OCO 122/20']},
                        ]

    ref_bond_coeffs = [
        # ua_dodecane
        {1: [130.0, 1.526, ' 1  SCP-SCS 6-2'],
         2: [130.0, 1.526, ' 2  SCS-SCS 2-2']},
        # ketene
        {1: [350.0, 1.305, ' 1  CHK-COK 47-110'],
         2: [170.0, 1.08, ' 2  CHK-HCK 46-47'],
         3: [350.0, 1.171, ' 3  COK-OCK 4-110']},
        # oxirene
        {1: [274.5, 1.34, ' 1  COO-COO 47-47'],
         2: [170.0, 1.08, ' 2  COO-HCO 46-47'],
         3: [225.0, 1.37, ' 3  COO-OCO 20-47']},
        # ethynol
        {1: [575.0, 1.21, ' 1  CHE-COE 19-19'],
         2: [210.0, 1.08, ' 2  CHE-HCE 19-46'],
         3: [275.0, 0.137, ' 3  COE-OHE 5-19 5-47'],
         4: [276.5, 0.945, ' 4  HOE-OHE 5-7']},
        # oxirene_bare
        {1: [' 1  COO-COO 47-47'],
         2: [' 2  COO-HCO 46-47'],
         3: [' 3  COO-OCO 20-47']},
                       ]

    ref_angle_coeffs = [
        # ua_dodecane
        {1: [31.5, 112.4,
             ' 1  SCP-SCS-SCS 6-2-2'],
         2: [31.5, 112.4,
             ' 2  SCS-SCS-SCS 2-2-2']},
        # ketene
        {1: [80.0, 180.0,
             ' 1  CHK-COK-OCK 4-110-47'],
         2: [17.5, 117.0,
             ' 2  COK-CHK-HCK 46-47-46'],
         3: [17.5, 117.0,
             ' 3  HCK-CHK-HCK 46-47-46']},
        # oxirene
        {1: [17.5, 120.0,
             ' 1  COO-COO-HCO 46-47-47'],
         2: [35.0, 123.0,
             ' 2  COO-COO-OCO 20-47-47'],
         3: [37.5, 60.0,
             ' 3  COO-OCO-COO 47-20-47 modified to smaller angle'],
         4: [17.5, 114.5,
             ' 4  HCO-COO-OCO 20-47-46']},
        # ethynol
        {1: [75.0, 180.0,
             ' 1  CHE-COE-OHE 5-19-19 13-19-19'],
         2: [56.0, 180.0,
             ' 2  COE-CHE-HCE 19-19-46'],
         3: [17.5, 109.0,
             ' 3  COE-OHE-HOE 7-5-19 7-5-47']},
        # oxirene_bare
        {1: [' 1  COO-COO-HCO 46-47-47'],
         2: [' 2  COO-COO-OCO 20-47-47'],
         3: [' 3  COO-OCO-COO 47-20-47 modified to smaller angle'],
         4: [' 4  HCO-COO-OCO 20-47-46']},
                        ]

    ref_dihedral_coeffs = [
        # ua_dodecane
        {1: [-3.4, 1.25, 3.1, 0.0,
             ' 1  SCP-SCS-SCS-SCS 6-2-2-2'],
         2: [-3.4, 1.25, 3.1, 0.0,
             ' 2  SCS-SCS-SCS-SCS 2-2-2-2']},
        # ketene
        {1: [0.0, 0.0, 0.0, 0.0,
             ' 1  HCK-CHK-COK-OCK 46-47-110-4']},
        # oxirene
        {1: [-3.5, 3.0, 0.0, 0.0,
             ' 1  COO-COO-OCO-COO 47-20-47-47 13-20-47-47'],
         2: [0.0, 14.0, 0.0, 0.0,
             ' 2  HCO-COO-COO-HCO 46-47-47-46'],
         3: [0.0, 14.0, 0.0, 0.0,
             ' 3  HCO-COO-COO-OCO 20-47-47-46'],
         4: [0.0, 0.0, 0.76, 0.0,
             ' 4  HCO-COO-OCO-COO 47-20-47-46 13-20-47-36'],
         5: [0.0, 14.0, 0.0, 0.0,
             ' 5  OCO-COO-COO-OCO 20-47-47-20 0-47-47-0']},
        # ethynol
        {1: [0.0, 0.0, 0.0, 0.0,
             ' 1  CHE-COE-OHE-HOE 7-5-19-19 0-19-19-0'],
         2: [0.0, 0.0, 0.0, 0.0,
             ' 2  HCE-CHE-COE-OHE 5-19-19-46 0-19-19-0']},
        # oxirene_bare
        {1: [' 1  COO-COO-OCO-COO 47-20-47-47 13-20-47-47'],
         2: [' 2  HCO-COO-COO-HCO 46-47-47-46'],
         3: [' 3  HCO-COO-COO-OCO 20-47-47-46'],
         4: [' 4  HCO-COO-OCO-COO 47-20-47-46 13-20-47-36'],
         5: [' 5  OCO-COO-COO-OCO 20-47-47-20 0-47-47-0']},
                           ]

    ref_improper_coeffs = [
        # ua_dodecane
        {},
        # ketene
        {1: [15.0, 180.0, ' 1  COK-HCK-CHK-HCK 110-46-47-46']},
        # oxirene
        {1: [15.0, 180.0, ' 1  COO-HCO-COO-OCO 47-46-47-20']},
        # ethynol
        {},
        # oxirene_bare
        {1: [' 1  COO-HCO-COO-OCO 47-46-47-20']},
                           ]

    @pytest.mark.parametrize("datafile, types",
                             zipRefs(ref_data_files_large(),
                                     ref_atom_types))
    def test_getAtomType(self, datafile, types):
        assert atoms.getAtomType(datafile) == types

    @pytest.mark.parametrize("datafile, atomdata, masses",
                             zipRefs(ref_data_files_large(),
                                     ref_atom_data,
                                     ref_masses))
    def test_getAtomData(self, datafile, atomdata, masses):
        assert atoms.getAtomData(datafile)[0][1] == atomdata
        assert atoms.getAtomData(datafile)[1] == masses

    @pytest.mark.parametrize("datafile, \
                             atomdata, \
                             masses, \
                             boxsizes, \
                             paircoeffs, \
                             bondcoeffs, \
                             anglecoeffs, \
                             dihedralcoeffs, \
                             impropercoeffs",
                             zipRefs(ref_data_files_large(),
                                     ref_atom_data,
                                     ref_masses,
                                     ref_boxsizes,
                                     ref_pair_coeffs,
                                     ref_bond_coeffs,
                                     ref_angle_coeffs,
                                     ref_dihedral_coeffs,
                                     ref_improper_coeffs))
    def test_getAllAtomDataLarge(self,
                                 datafile,
                                 atomdata,
                                 masses,
                                 boxsizes,
                                 paircoeffs,
                                 bondcoeffs,
                                 anglecoeffs,
                                 dihedralcoeffs,
                                 impropercoeffs):
        assert atoms.getAllAtomData(datafile)[0][1] == atomdata
        assert atoms.getAllAtomData(datafile)[5] == masses
        assert atoms.getAllAtomData(datafile)[6] == boxsizes
        assert atoms.getAllAtomData(datafile)[8] == paircoeffs
        assert atoms.getAllAtomData(datafile)[9] == bondcoeffs
        assert atoms.getAllAtomData(datafile)[10] == anglecoeffs
        assert atoms.getAllAtomData(datafile)[11] == dihedralcoeffs
        assert atoms.getAllAtomData(datafile)[12] == impropercoeffs


class TestAtomDataSmall:
    ref_bonds = [
                 # ketene
                 {1: [1, 1, 2],
                  2: [2, 1, 4],
                  3: [2, 1, 5],
                  4: [3, 2, 3]},
                 # oxirene
                 {1: [1, 1, 2],
                  2: [3, 1, 3],
                  3: [2, 1, 4],
                  4: [3, 2, 3],
                  5: [2, 2, 5]},
                 # ethynol
                 {1: [2, 1, 4],
                  2: [1, 1, 2],
                  3: [3, 2, 3],
                  4: [4, 3, 5]},
                 # oxirene_bare
                 {1: [1, 1, 2],
                  2: [3, 1, 3],
                  3: [2, 1, 4],
                  4: [3, 2, 3],
                  5: [2, 2, 5]},
                 ]

    ref_angles = [
        # ketene
        {1: [2, 2, 1, 4],
         2: [2, 2, 1, 5],
         3: [3, 4, 1, 5],
         4: [1, 1, 2, 3]},
        # oxirene
        {1: [2, 2, 1, 3],
         2: [1, 2, 1, 4],
         3: [4, 3, 1, 4],
         4: [2, 1, 2, 3],
         5: [1, 1, 2, 5],
         6: [4, 3, 2, 5],
         7: [3, 1, 3, 2]},
        # ethynol
        {1: [2, 2, 1, 4],
         2: [1, 1, 2, 3],
         3: [3, 2, 3, 5]},
        # oxirene_bare
        {1: [2, 2, 1, 3],
         2: [1, 2, 1, 4],
         3: [4, 3, 1, 4],
         4: [2, 1, 2, 3],
         5: [1, 1, 2, 5],
         6: [4, 3, 2, 5],
         7: [3, 1, 3, 2]},
                  ]

    ref_dihedrals = [
        # ketene
        {1: [1, 4, 1, 2, 3],
         2: [1, 5, 1, 2, 3]},
        # oxirene
        {1: [5, 3, 1, 2, 3],
         2: [3, 3, 1, 2, 5],
         3: [3, 4, 1, 2, 3],
         4: [2, 4, 1, 2, 5],
         5: [1, 2, 1, 3, 2],
         6: [4, 4, 1, 3, 2],
         7: [1, 1, 2, 3, 1],
         8: [4, 5, 2, 3, 1]},
        # ethynol
        {1: [2, 4, 1, 2, 3],
         2: [1, 1, 2, 3, 5]},
        # oxirene_bare
        {1: [5, 3, 1, 2, 3],
         2: [3, 3, 1, 2, 5],
         3: [3, 4, 1, 2, 3],
         4: [2, 4, 1, 2, 5],
         5: [1, 2, 1, 3, 2],
         6: [4, 4, 1, 3, 2],
         7: [1, 1, 2, 3, 1],
         8: [4, 5, 2, 3, 1]},
                     ]

    ref_impropers = [
        # ketene
        {1: [1, 5, 1, 4, 2]},
        # oxirene
        {1: [1, 3, 1, 4, 2],
         2: [1, 3, 2, 5, 1]},
        # ethynol
        {},
        # oxirene_bare
        {1: [1, 3, 1, 4, 2],
         2: [1, 3, 2, 5, 1]},
                     ]

    @pytest.mark.parametrize("datafile, \
                             bonds, \
                             angles, \
                             dihedrals, \
                             impropers",
                             zipRefs(ref_data_files(),
                                     ref_bonds,
                                     ref_angles,
                                     ref_dihedrals,
                                     ref_impropers))
    def test_getAllAtomData2(self,
                             datafile,
                             bonds,
                             angles,
                             dihedrals,
                             impropers):
        assert atoms.getAllAtomData(datafile)[1] == bonds
        assert atoms.getAllAtomData(datafile)[2] == angles
        assert atoms.getAllAtomData(datafile)[3] == dihedrals
        assert atoms.getAllAtomData(datafile)[4] == impropers


class TestRanges:
    """
    first, last
    first, ~mid
    ~mid, last
    ~third, ~third
    """

    ref_AtomRanges = [
        # ua_dodecane
        [2400, 1, 3, 2398, 2396, 1, 0, 0],
        # ketene
        [5, 1, 3, 3, 1, 1, 0, 0],
        # oxirene
        [5, 1, 3, 3, 1, 1, 0, 0],
        # ethynol
        [5, 1, 3, 3, 1, 1, 0, 0],
        # oxirene_bare
        [5, 1, 3, 3, 1, 1, 0, 0],
        ]

    @pytest.mark.parametrize("datafile, result",
                             zipRefs(ref_data_files_large(), ref_AtomRanges))
    @pytest.mark.parametrize(
        "idx, start, stop", [
            (0, 'first', 'last'),
            (1, 'first', 1),
            (2, 'first', 3),
            (3, 3, 'last'),
            (4, 5, 'last'),
            (5, 2, 2),
            (6, 3, 2),
            (7, 3E10, 3),
        ])
    def test_getAtomRange(self, datafile, result, idx, start, stop):
        atomdata = atoms.getAllAtomData(datafile)[0]
        assert len(atoms.getAtomRange(atomdata, start, stop)) == result[idx]
