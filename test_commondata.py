import pytest
import commondata_p3 as cdp3


@pytest.mark.parametrize(
    "filename, lines", [
        ("Fe2O3_50_down.lammps", 14120),
        ("Fe2O3_50_up.lammps", 14120),
        ("Fe2O3_100_down.lammps", 56360),
        ("Fe2O3_100_up.lammps", 56360)
    ])
def test_readAll(filename, lines):
    assert len(cdp3.readAll(filename)) == lines


@pytest.mark.parametrize(
    "filename, lines", [
        ("Fe2O3_50_down.lammps.gz", 14120),
        ("Fe2O3_50_up.lammps.gz", 14120),
        ("Fe2O3_100_down.lammps.gz", 56360),
        ("Fe2O3_100_up.lammps.gz", 56360)
    ])
def test_readAllGzip(filename, lines):
    assert len(cdp3.readAllGzip(filename)) == lines


@pytest.mark.parametrize(
    "lines, result",
    [(
        ["not this line",
         "10000",
         "ITEM: BOXLENGTH",
         "100",
         "200",
         "300",
         "ITEM: NUMBER OF ATOMS",
         "16889"],
        16889)
     ])
def test_getNatoms(lines, result):
    assert cdp3.getNatoms(lines) == result
