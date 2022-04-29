import pytest
import lmptools.mergedatafiles_p3 as mdf3


@pytest.mark.parametrize("testargs, message", [
    (["mergedatafile_p3.py"], 'Missing settings'),
    (["mergedatafile_p3.py", "merge2.yaml", "mistake"], 'This script'),
    ])
def test_readInputFile_fail(testargs, message):
    try:
        with pytest.raises(mdf3.MissingSettingsFile) as excinfo:
            mdf3.readInputFile(testargs)
            assert message in str(excinfo.value)
    except SystemExit:
        pass


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
