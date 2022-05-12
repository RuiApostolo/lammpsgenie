# lmptools

This package is a collection of tools for handling LAMMPS data and dump files.

## Script

The main way to interact with lmptools is using the script `mergedatafiles`, which should be called with a settings file (yaml format) like so:

```
mergedatafiles merge_settings.yaml
```

If called without any argument, it will search the current directory for a `merge.yaml` and `merge.yml` (in this order) and use the first one it finds, or return an error.


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




[numpy styleguide]: https://numpydoc.readthedocs.io/en/latest/format.html
