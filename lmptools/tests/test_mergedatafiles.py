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


@pytest.mark.parametrize("testargs, message", [
    (["mergedatafile_p3.py"], 'Missing settings'),
    (["mergedatafile_p3.py", "merge2.yaml", "mistake"], 'This script'),
    ])
def test_readInputFile_fail(testargs, message):
    try:
        with pytest.raises(mdf3.MissingSettingsFile) as excinfo:
            mdf3.readInputFile(testargs)
            assert message in str(excinfo.value)
    # catch systemexit and prevent from failing test
    except SystemExit:
        pass


class TestSettings:
    ref_mergeYaml = [
        {"filename": "Fe2O3_50_down.lammps",
         "minormax": "max",
         "value": -1.5},
        {"filename": "Fe2O3_50_up.lammps",
         "minormax": "min",
         "value": 80.0},
        {"filename": "uadodecane.data",
         "minormax": "min",
         "value": 1.5},
        ]

    ref_minmax = {
        # each line: x_min, y_min, z_min, x_max, y_max, z_max
        "Fe2O3_50_down.lammps": {
            min: {'x': 0.303, 'y': 0.0, 'z': -8.606},
            max: {'x': 54.785, 'y': 49.609001, 'z': 0.0},
            },
        "Fe2O3_50_up.lammps": {
            min: {'x': 0.303, 'y': 0.0, 'z': 458},
            max: {'x': 54.785, 'y': 49.609001, 'z': 466.606},
            },
        "uadodecane.data": {
            min: {'x': 1.00002, 'y': 1.000255, 'z': 0.943834},
            max: {'x': 44.061073, 'y': 44.005527, 'z': 43.997902},
            },
    }

    def test__openYaml(self, base_path, monkeypatch):
        monkeypatch.chdir("tests")
        assert mdf3._openYaml("merge.yaml") == self.ref_mergeYaml

    @pytest.mark.parametrize("args", [
        (["", "merge.yaml"]),
        # testing default value:
        ([""]),
        ])
    def test_readInputFile_success(self,
                                   base_path,
                                   monkeypatch,
                                   args):
        monkeypatch.chdir("tests")
        settings = mdf3.readInputFile(args)
        assert settings == self.ref_mergeYaml

    @pytest.mark.parametrize("file", zipRefs(ref_mergeYaml))
    @pytest.mark.parametrize("coordinate", ['x', 'y', 'z'])
    @pytest.mark.parametrize("function", [min, max])
    def test__limit(self,
                    base_path,
                    monkeypatch,
                    file,
                    function,
                    coordinate):
        monkeypatch.chdir("tests")
        topology = mdf3.readTopology(file[0])
        assert mdf3._limit(function, coordinate, topology) == pytest.approx(
            self.ref_minmax[file[0]['filename']][function][coordinate])

    def test_readTopologies(self, base_path, monkeypatch):
        monkeypatch.chdir("tests")
        topologies = mdf3.readTopologies(self.ref_mergeYaml)
        assert len(topologies) == 3

#  class TestTopologies:
