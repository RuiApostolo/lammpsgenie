from yaml import full_load, dump
import lmptools.mergedatafiles_p3 as mdf3


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
    mdf3.main(args)
    with open('expected_merged_uadodecane.lammps', 'r') as f:
        expected = f.readlines()
    with open(out_file, 'r') as f:
        result = f.readlines()
    assert result[1:] == expected[1:]
