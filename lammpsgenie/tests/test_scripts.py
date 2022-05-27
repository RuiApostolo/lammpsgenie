from yaml import full_load, dump
from pathlib import Path
import lammpsgenie.mergedatafiles as mdf
import lammpsgenie.copy_ironfiles as cif


def modified_settings(orig_file, new_path):
    with open(orig_file, 'r') as f:
        settings = full_load(f)
    settings['output']['filename'] = str(new_path)
    return settings


def test_mergedatafiles_script(mock_path, tmp_path):
    # name setup
    settings_name = tmp_path / 'merge.yaml'
    out_file = tmp_path / 'merged_uadodecane.lammps'
    with open(settings_name, 'w') as settings_file:
        dump(modified_settings('merge.yaml', out_file), settings_file)

    args = ['mergedatafiles_p3.py', str(settings_name)]
    mdf.main(args)
    with open('expected_merged_uadodecane.lammps', 'r') as f:
        expected = f.readlines()
    with open(out_file, 'r') as f:
        result = f.readlines()
    assert result[1:] == expected[1:]


def test_saveiron_script(tmp_path, monkeypatch):
    iron_50 = ['Fe2O3_50_down.lammps', 'Fe2O3_50_up.lammps']
    iron_100 = ['Fe2O3_100_down.lammps', 'Fe2O3_100_up.lammps']
    iron_all = iron_50 + iron_100

    # hack. there's likely a much better way to solve this
    cwd = Path(__file__).parent.resolve()
    monkeypatch.chdir(tmp_path)
    cif.save_iron_all()
    monkeypatch.chdir(str(cwd) + "/../data")
    for file in iron_all:
        with open(file, 'r') as f:
            expected = f.readlines()
        with open(tmp_path/file, 'r') as f:
            result = f.readlines()
        assert result == expected
