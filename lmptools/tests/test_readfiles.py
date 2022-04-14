import pytest
from conftest import zipRefs, ref_data_files_large
import lmptools.readfiles as rdfl


class TestReadAll:
    ref_number_lines = [8452, 80, 94, 75]

    @pytest.mark.parametrize("filename, lines",
                             zipRefs(ref_data_files_large(),
                                     ref_number_lines))
    def test_readAll_fromParameters(self, filename, lines):
        result = len(rdfl.readAll(filename))
        assert result == lines

    @pytest.mark.parametrize(
        "filename, nlines", [
            ("tests/dump.uadodecane.lammpstrj", 14454)
        ])
    def test_readAll_fromFile(self, filename, nlines):
        assert len(rdfl.readAll(filename)) == nlines

    @pytest.mark.parametrize("filename, lines",
                             zipRefs(ref_data_files_large(),
                                     ref_number_lines))
    def test_readAllGzip(self, filename, lines):
        assert len(rdfl.readAllGzip(filename+".gz")) == lines
