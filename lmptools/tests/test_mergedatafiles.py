import pytest
import lmptools.mergedatafiles_p3 as mdf3
from conftest import zipRefs
from copy import deepcopy


ref_mergeYamlOut = {
    "filename": "merged_uadodecane.lammps",
    "boxsize": {
        "xlo": 0.0,
        "xhi": 55.088,
        "ylo": 0.0,
        "yhi": 50.38,
        "zlo": -20.0,
        "zhi": 120.0,
    }
}

ref_mergeYamlIn = {
    "Fe2O3_50_down.lammps": {
        "minormax": max,
        "value": -1.5},
    "Fe2O3_50_up.lammps": {
        "minormax": min,
        "value": 80.0},
    "uadodecane.data": {
        "minormax": min,
        "value": 1.5},
    }

ref_mergeYamlAll = {
    **{'outputs': ref_mergeYamlOut},
    **{'inputs': ref_mergeYamlIn}
    }


@pytest.mark.parametrize("test_dict, message", [
    ({"name": {"minormax": "man", "value": 30.3}}, "man"),
    ({"name": {"minormax": "", "value": 30.3}}, ""),
    ])
def test__strToFunction(test_dict, message):
    with pytest.raises(ValueError) as excinfo:
        mdf3._strToFunction(test_dict)
    assert message in str(excinfo.value)


@pytest.mark.parametrize("start, expected, count", [
    ({1: [1.002, -0.3, ' SCP'], 2: [2.0, 0.0, 'SCS']},
     {6: [1.002, -0.3, ' SCP'], 7: [2.0, 0.0, 'SCS']},
     5),
    ({1: {'mol': 1,
          'type': 'SCP',
          'charge': -0.5,
          'x': 1.0,
          'y': 2.0,
          'z': 2.0}},
     {2: {'mol': 1,
          'type': 'SCP',
          'charge': -0.5,
          'x': 1.0,
          'y': 2.0,
          'z': 2.0}},
     1),
])
def test__shiftKey(start, expected, count):
    assert mdf3._shiftKey(start, count) == expected


@pytest.mark.parametrize("dict_a, dict_b, expected", [
    ({'atoms': 1},
     {'atoms': 2},
     {'atoms': 3}),
    ({'atoms': 10, 'bonds': 20, 'angles': 30},
     {'atoms': 21, 'bonds': 33, 'angles': 47},
     {'atoms': 31, 'bonds': 53, 'angles': 77}),
    ({'atoms': 13, 'bonds': 45, 'angles': 35},
     {'atoms': 21, 'bonds': 26, 'dihedrals': 46},
     {'atoms': 34, 'bonds': 71, 'angles': 35, 'dihedrals': 46}),
])
def test__combineDicts(dict_a, dict_b, expected):
    assert mdf3._combineDicts(dict_a, dict_b) == expected


@pytest.mark.parametrize("dict_a, dict_b, expected", [
    ({'atoms': 10, 'bonds': 20, 'angles': 30},
     {'angles': 30, 'dihedrals': 40},
     {'atoms': 10, 'bonds': 20, 'angles': 30, 'dihedrals': 40}),
])
def test__mergeDictsPass(dict_a, dict_b, expected):
    assert mdf3._mergeDicts(dict_a, dict_b) == expected


@pytest.mark.parametrize("dict_a, dict_b, message", [
    ({'atoms': 10, 'bonds': 20, 'angles': 30},
     {'angles': 30, 'bonds': 90, 'dihedrals': 40},
     'bonds')
])
def test__mergeDictsRaise(dict_a, dict_b, message):
    with pytest.raises(mdf3.ValueExists) as excinfo:
        mdf3._mergeDicts(dict_a, dict_b)
    assert message in str(excinfo.value)


@pytest.mark.parametrize("testargs, message", [
    (["mergedatafile_p3.py"], 'Missing settings'),
    (["mergedatafile_p3.py", "merge2.yaml", "mistake"], 'This script'),
    ])
def test_readInputFile_fail(testargs, message):
    with pytest.raises(mdf3.MissingSettingsFile) as excinfo:
        mdf3.readInputFile(testargs)
    assert message in str(excinfo.value)


class TestSettings:

    ref_minmax = {
        # each line: x_min, y_min, z_min, x_max, y_max, z_max
        "Fe2O3_50_down.lammps": {
            min: {'x': 0.303, 'y': 0.0, 'z': -8.606},
            max: {'x': 54.785, 'y': 49.609001, 'z': 0.0},
            },
        "Fe2O3_50_up.lammps": {
            min: {'x': 0.303, 'y': 0.0, 'z': 458.0},
            max: {'x': 54.785, 'y': 49.609001, 'z': 466.606},
            },
        "uadodecane.data": {
            min: {'x': 1.00002, 'y': 1.000255, 'z': 0.943834},
            max: {'x': 44.061073, 'y': 44.005527, 'z': 43.997902},
            },
    }

    def approx_nested_dict(self, dictionary):
        # uses pytest.approx() to values in nested dicts
        result = {}
        for key in dictionary:
            result[key] = pytest.approx(dictionary[key])
        return result

    def approx_double_nested_dict(self, dictionary):
        # uses pytest.approx() to values in double-nested dicts
        result = {}
        for key in dictionary:
            result[key] = self.approx_nested_dict(dictionary[key])
        return result

    def test__openYaml(self, mock_path):
        assert mdf3._openYaml("merge.yaml") == \
            (ref_mergeYamlIn, ref_mergeYamlOut)

    @pytest.mark.parametrize("args", [
        (["", "merge.yaml"]),
        # testing default value:
        ([""]),
        ])
    def test_readInputFile_success(self,
                                   mock_path,
                                   args):
        settings, outputfile = mdf3.readInputFile(args)
        assert settings == ref_mergeYamlIn
        assert outputfile == ref_mergeYamlOut

    @pytest.mark.parametrize("file", zipRefs(ref_mergeYamlIn))
    @pytest.mark.parametrize("coordinate", ['x', 'y', 'z'])
    @pytest.mark.parametrize("function", [min, max])
    def test__limit(self,
                    mock_path,
                    file,
                    function,
                    coordinate):
        topology = mdf3.readTopology(file[0])
        assert mdf3._limit(function, coordinate, topology) == pytest.approx(
            self.ref_minmax[file[0]][function][coordinate])


class TestTopologies(TestSettings):
    @pytest.fixture
    def topologies(self):
        return mdf3.readTopologies(ref_mergeYamlIn)

    ref_limits = {min: {
                    'x': 0.303,
                    'y': 0.0,
                    'z': -8.606},
                  max: {
                    'x': 54.785,
                    'y': 49.609001,
                    'z': 466.606}}

    ref_newminmax = {
        # each line: x_min, y_min, z_min, x_max, y_max, z_max
        "Fe2O3_50_down.lammps": {
            min: {'x': 0.303, 'y': 0.0, 'z': -10.106},
            max: {'x': 54.785, 'y': 49.609001, 'z': -1.5},
            },
        "Fe2O3_50_up.lammps": {
            min: {'x': 0.303, 'y': 0.0, 'z': 80.0},
            max: {'x': 54.785, 'y': 49.609001, 'z': 88.606},
            },
        "uadodecane.data": {
            min: {'x': 1.00002, 'y': 1.000255, 'z': 1.5},
            max: {'x': 44.061073, 'y': 44.005527, 'z':  44.554068},
            },
    }

    def test_readTopologies(self, mock_path, topologies):
        assert len(topologies) == 3

    def test_limitsAllTopologies(self, mock_path, topologies):
        assert mdf3.limitsAllTopologies(topologies) == \
            self.approx_double_nested_dict(self.ref_minmax)

    def test_AbsoluteLimitsTopologies(self, mock_path, topologies):
        assert mdf3.absoluteLimitsTopologies(topologies) == \
            self.approx_nested_dict(self.ref_limits)


class TestMerge(TestTopologies):
    @pytest.fixture
    def shifted_topologies(self, topologies):
        return mdf3.shiftTopologies(topologies,
                                    self.ref_minmax,
                                    ref_mergeYamlIn)

    @pytest.fixture
    def expected_merged(self):
        return mdf3.readTopology('expected_merged_uadodecane.lammps')

    @pytest.fixture
    def merged_undertest1(self, shifted_topologies):
        return mdf3.mergeTopologies(
            shifted_topologies, ref_mergeYamlOut['boxsize'])

    @pytest.fixture
    def merged_undertest2(self, merged_undertest1):
        merged_undertest2 = deepcopy(merged_undertest1)
        merged_undertest2['masses'] = {}
        merged_undertest2['atomdata'] = {}
        return merged_undertest2

    def test_shiftTopologies(self, mock_path, topologies, shifted_topologies):
        assert topologies != shifted_topologies

    def test_limitsAllTopologies(self,
                                 mock_path,
                                 topologies,
                                 shifted_topologies):
        assert mdf3.limitsAllTopologies(shifted_topologies) == \
            self.approx_double_nested_dict(self.ref_newminmax)

    def test_mergeTopologies(self,
                             mock_path,
                             shifted_topologies,
                             merged_undertest1,
                             expected_merged):
        for atom in merged_undertest1['atomdata']:
            assert merged_undertest1['atomdata'][atom] == \
                pytest.approx(expected_merged['atomdata'][atom])
        for propert in [prop for prop in merged_undertest1
                        if prop != 'atomdata']:
            assert merged_undertest1[propert] == expected_merged[propert]

    @pytest.mark.parametrize("topology, result", [
        (pytest.lazy_fixture('merged_undertest1'),
         'expected_merged_uadodecane.lammps'),
        (pytest.lazy_fixture('merged_undertest2'),
         'expected_merged_cut_uadodecane.lammps'),
    ])
    def test_writeTopology(self, mock_path, tmpdir, topology, result):
        file = tmpdir.join('out.lammps')
        mdf3.writeTopology(topology, file)
        with open(result, 'r') as f:
            expected = f.readlines()
        assert file.readlines()[1:] == expected[1:]


@pytest.mark.parametrize("dictionary, keys, expected", [
    ({1: {1: [1, 2, 3, 4], 2: [3, 4, 5, 6]},
      2: {1: [10, 11, 12], 2: [11, 15, 17]},
      },
     {1: 'One', 2: 'Two'},
     '\nOne\n\n1 1 2 3 4\n2 3 4 5 6\n\nTwo\n\n1 10 11 12\n2 11 15 17\n',
     ),
    ({1: {1: [1, 2, 3, 4], 2: [3, 4, 5, 6]},
      2: {}
      },
     {1: 'One', 2: 'Two'},
     '\nOne\n\n1 1 2 3 4\n2 3 4 5 6\n',
     ),
])
def test__getStringBADI(dictionary, keys, expected):
    assert mdf3._getStringBADI(dictionary, keys) == expected
