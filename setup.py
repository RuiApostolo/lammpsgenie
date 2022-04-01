import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lmptools",
    version="0.0.2",
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
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6"
)
