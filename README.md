# lmptools

This package is a collection of tools for handling LAMMPS data and dump files.

## Scripts

### mergedatafiles

The main way to interact with lmptools is using the script `mergedatafiles`, which should be called with a settings file (YAML format) like so:

```
mergedatafiles merge_settings.yaml
```

If called without any argument, it will search the current directory for a `merge.yaml` and `merge.yml` (in this order) and use the first one it finds, or return an error.

The YAML file should contain two dictionaries (`inputs` and `output`), with the following format:

```
output:
  filename: <filename/path where to write merged file>
  boxsize: # new boxsize for the merged file
    xlo: float
    xhi: float
    ylo: float
    yhi: float
    zlo: float
    zhi: float
inputs:
  # as many of the following as needed
  <filename/path of a data file to read from>:
    minormax: <min or max>
    value: float
```

The script will gather the information from all the input files, shift coordinates to make sure they match the min/max values as given, and then write a new merged data file, with the new box limits provided.

Example:

```
output:
  filename: merged_uadodecane.lammps
  boxsize:
    xlo: 0.0
    xhi: 55.088
    ylo: 0.0
    yhi: 50.38
    zlo: -20.0
    zhi: 120.0
inputs:
  ../data/Fe2O3_50_down.lammps:
    minormax: max
    value: -1.5
  ../data/Fe2O3_50_up.lammps:
    minormax: min
    value: 80.0
  uadodecane.data:
    minormax: min
    value: 1.5
```

### saveiron

There are three variants of the `saveiron` script that write to the current working directory data files of the iron oxide surface:

* saveiron50
* saveiron100
* saveironall

The first one writes the files of surfaces with side 50Å, the second the files for surfaces with side 100Å, and the last one writes the files for surfaces of both 50 and 100Å.


## Development

### Virtual Environment

A virtual environment (venv) is like a sandbox, lets you install fresh version of packages (including local ones) without conflicting with your main python installation.
It is recommended to use `virtualenvwrapper`, a series of shell scripts that automates the creation, activation, deactivation, (among other things) of venvs.
To install `virtualenvwrapper` use the command:

```
pip install virtualenvwrapper
```

Then you might need to export some variables on your `~/.profile` (order matters):

```
# virtualenvwrapper
# if you have python2 installed, you might need:
export VIRTUALENVWRAPPER_PYTHON=$(which python3)
# location of virtualenvwrapper.sh
source $(which virtualenvwrapper.sh)
# where to store the files for each venv
export WORKON_HOME=$HOME/.virtualenvs
# used when creating new projects
export PROJECT_HOME=$HOME/code
```

And restart your shell (or source `~/.profile`).

To create a new venv, use:

```
mkvirtuaenv <name>
```

This will create the venv and activate it.
Exit the venv with

```
deactivate
```

and return to the venv with

```
workon <name>
```

Sometimes, you need to add the path of the project to a venv.
That is done, when one has the venv activated, with:

```
add2virtualenv /absolute/path/to/package/
```

### Tests

You will need pytest and some plugins to run tests.
The required python packages are provided in test\_requirements.txt, and can be installed with

```
pip install -r test_requirements.txt
```

This is better done inside a virtual environment, see the previous section.

Tests can be run from the package folder with:

```
pytest
```

This will run every test in the `tests/` folder.
To run only one of the test files, use:

```
pytest tests/test_<name>.py
```

During development, a good way to keep an eye on tests when changing code (or adding new tests/test data) is to keep a separate terminal running `pytest-watch`.
This will scan the module and test files, and run the new test, or the tests affected by any changes, once it detects a file has been written.
Use the commmand:

```
ptw --runner "pytest --testmon"
```

To check that the test cover the entirety of your code, use the coverage plugin for pytest:

```
pytest --cov=lmptools/
```

There are several settings that can be changed for this tool, to do so, add a `.coveragerc` file to the project directory.
The following are a good set of defaults:

```
[run]
branch = True

[paths]
source =
  lmptools

[report]
show_missing = True
skip_empty = True
precision = 2
```

### DocStrings

The [numpy styleguide] was followed for DocString documentation.


## Release

Install required packages with:

```
pip install --upgrade build twine
```

1. Modify code
2. Assure tests pass, and cover 100% of the code.
3. Modify version in `_version.py`
4. Create a build with `python3 -m build`.
5. Check that the distribution files pass checks with `twine check dist/*`
6. Upload to PyPi with `python3 -m twine upload --repository pypi dist/lmptools-<version>*`


Note: you need an account on pypi, and the necessary rights to upload, and a [registered token] saved on `.pypirc`

To test your distribution, you might want to test upload to the PyPi test repository with `python3 -m twine upload --repository pypi dist/lmptools-<version>*`
(Needs a separate registered account).

[numpy styleguide]: https://numpydoc.readthedocs.io/en/latest/format.html
[registered token]: https://pypi.org/help/#apitoken
