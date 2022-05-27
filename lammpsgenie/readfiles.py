import gzip
"""
Functions to read from files
"""

__all__ = [
    'readAll',
    'readAllGzip',
]


def readAll(filename):
    """
    Reads the entire file in line-by-line and returns list of lines.

    Parameters
    ----------
    filename : str
        Name of file to be read.

    Returns
    -------
    lines : list
        A list containing the lines of the given file.
    """

    with open(filename, "r") as ifile:
        lines = ifile.read().splitlines()
    return lines


def readAllGzip(filename):
    """
    Reads the entire gzipped file in line-by-line and returns list of
    lines.

    Parameters
    ----------
    filename : str
        Name of file to be read.

    Returns
    -------
    lines: list of str
        A list containing the lines of given file.
    """

    with gzip.open(filename, "r") as ifile:
        lines = ifile.readlines()
    return lines
