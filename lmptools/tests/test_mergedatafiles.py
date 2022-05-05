import pytest
import lmptools.mergedatafiles_p3 as mdf3
from conftest import zipRefs


@pytest.mark.parametrize("message, code", [
    ("Message 1", 3),
    ("123 Message test", 4),
    ])
def test__myExit(message, code, capsys):
    with pytest.raises(SystemExit) as exc_info:
        mdf3._myExit(message, code)
    assert exc_info.type is SystemExit
    assert exc_info.value.code == code
    captured = capsys.readouterr()
    assert captured.out == message + "\n"


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
def test__shiftBy(start, expected, count):
    assert mdf3._shiftBy(start, count) == expected


@pytest.mark.parametrize("dict_a, dict_b, expected", [
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
    @pytest.fixture
    def mock_path(self, monkeypatch):
        monkeypatch.chdir("tests")

    ref_mergeYamlOut = "merged.data.lammps"

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
            (self.ref_mergeYamlIn, self.ref_mergeYamlOut)

    @pytest.mark.parametrize("args", [
        (["", "merge.yaml"]),
        # testing default value:
        ([""]),
        ])
    def test_readInputFile_success(self,
                                   mock_path,
                                   args):
        settings, outputfile = mdf3.readInputFile(args)
        assert settings == self.ref_mergeYamlIn
        assert outputfile == self.ref_mergeYamlOut

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
        return mdf3.readTopologies(self.ref_mergeYamlIn)

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

    def test_shiftTopologies(self, mock_path, topologies):
        shifted_topologies = mdf3.shiftTopologies(topologies,
                                                  self.ref_minmax,
                                                  self.ref_mergeYamlIn)
        assert topologies != shifted_topologies
        assert mdf3.limitsAllTopologies(shifted_topologies) == \
            self.approx_double_nested_dict(self.ref_newminmax)
