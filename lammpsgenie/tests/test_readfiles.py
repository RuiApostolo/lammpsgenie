import pytest
from conftest import zipRefs, ref_all_data_fs
import lmptools.readfiles as rdfl


class TestReadAll:
    ref_nlines_datafs = [8452, 80, 94, 75, 94]
    ref_nlines_ironfs = [14120, 14120, 56360, 56360]

    @pytest.mark.parametrize("filename, lines",
                             zipRefs(ref_all_data_fs(),
                                     ref_nlines_datafs + ref_nlines_ironfs))
    def test_readAll_Data(self, filename, lines):
        print(filename, lines)
        result = len(rdfl.readAll(filename))
        assert result == lines

    @pytest.mark.parametrize(
        "filename, nlines", [
            ("tests/dump.uadodecane.lammpstrj", 14454)
        ])
    def test_readAll_Dump(self, filename, nlines):
        assert len(rdfl.readAll(filename)) == nlines

    @pytest.mark.parametrize("filename, lines",
                             zipRefs(ref_all_data_fs(),
                                     ref_nlines_datafs + ref_nlines_ironfs))
    def test_readAllGzip(self, filename, lines):
        assert len(rdfl.readAllGzip(filename+".gz")) == lines
