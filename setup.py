import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lmptools",
    version="0.1.0",
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
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Development Status :: 3 - Alpha",
        "Operating System :: Unix"
        ],
    package_dir={"": "lmptools"},
    packages=setuptools.find_packages(where="lmptools"),
    python_requires=">=3.6",
    install_requires=[
        "PyYAML<=5.3.1"
    ],
    #  scripts=['lmptools/mergedatafiles_p3.py'],
    entry_points={
        'console_scripts': ['mergedatafiles=lmptools.mergedatafiles_p3:main'],
    }
)
