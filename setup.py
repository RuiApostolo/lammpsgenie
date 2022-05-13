from setuptools import setup, find_packages
from distutils.util import convert_path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


main_ns = {}
ver_path = convert_path("lmptools/_version.py")
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)


setup(
    name="lmptools",
    version=main_ns['__version__'],
    author="Rui Apóstolo",
    author_email="ruiapostolo@gmail.com",
    description="A package to help with LAMMPS data and dump files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RuiApostolo/lmpdtmrg",
    project_urls={
        "Bug Tracker": "https://github.com/RuiApostolo/lmpdtmrg/issues"
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        ],
    packages=find_packages(where="."),
    python_requires=">=3.6",
    install_requires=[
        "PyYAML>=5.3.1"
    ],
    entry_points={
        'console_scripts': [
            'mergedatafiles=lmptools.mergedatafiles:main',
            'saveiron50=lmptools.copy_ironfiles:save_iron_50',
            'saveiron100=lmptools.copy_ironfiles:save_iron_100',
            'saveironall=lmptools.copy_ironfiles:save_iron_all',
        ],
    },
    include_package_data=True,
    package_data={"": ["data/*.lammps"]},
)
