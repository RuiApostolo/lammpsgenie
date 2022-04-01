import pytest
import lmptools.readfiles as rdfl


@pytest.mark.parametrize(
    "filename, lines", [
        ("tests/Fe2O3_50_down.lammps", 14120),
        ("tests/Fe2O3_50_up.lammps", 14120),
        ("tests/Fe2O3_100_down.lammps", 56360),
        ("tests/Fe2O3_100_up.lammps", 56360)
    ])
def test_readAll_fromParameters(filename, lines, dummy_dump):
    result = len(rdfl.readAll(filename))
    assert result == lines


@pytest.mark.parametrize(
    "filename, nlines", [
        ("tests/dump.uadodecane.lammpstrj", 14454)
    ])
def test_readAll_fromFile(filename, nlines):
    assert len(rdfl.readAll(filename)) == nlines


@pytest.mark.parametrize(
    "filename, lines", [
        ("tests/Fe2O3_50_down.lammps.gz", 14120),
        ("tests/Fe2O3_50_up.lammps.gz", 14120),
        ("tests/Fe2O3_100_down.lammps.gz", 56360),
        ("tests/Fe2O3_100_up.lammps.gz", 56360)
    ])
def test_readAllGzip(filename, lines):
    assert len(rdfl.readAllGzip(filename)) == lines
