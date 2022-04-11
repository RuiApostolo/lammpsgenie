Test molecule generation and configuration
==========================================

Three molecules were generated, for testing purposes: ketene (1), oxirene (2), and ethynol (3).
They were chosen because they have the same number of atoms, but offer distinct geometries and number of atom types.
Values are not guaranteed to be accurate, but are sufficient for testing purposes.
Steps to generate all three molecules are shown below, with commands substituted with `<variables>` where appliable, and `<variable>` values are detailed elsewhere.

## Requirements

These instructions make use of third party tools:

1. OpenBabel: available in aptitude with `sudo apt install openbabel` or [Download](http://openbabel.org/wiki/Category:Installation)
2. VMD: [Download](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
3. VMD's packages: topotools and pbctools, which come by default on most VMD installations.

## Structure

1. Ketene - H2C=C=O
2. Oxirene - HC(-O-)=CH
3. Ethynol - HC#COH

## Smiles

1. C=C=O
2. C1=CO1
3. C#CO

To generate the 3d structure of a molecule and save it to a .xyz file, use the command:
```
obabel -i smi -:<smiles> -h --gen3d -o xyz -O <molecule>.xyz
```
where `<molecule>` is the name of the molecule.


To generate the lammps data file (using VMD), use the command:
```
vmd -dispdev text -e <molecule>.vmd
```
which contains the charge values for each atom, then set the pair, bond, angle, dihedral, and improper parameters:


## Force field parameters

All parameters used here are taken from OPLS-AA.


### Ketene

Structure with custom types: `HCK(2)-CHK=COK=OCK`

#### Atom Types

* HCK: 899/46
* CHK: 900/47
* COK: 904/110
* OCK: 905/4

#### Charges

* HCK: +0.1500
* CHK: -0.2500
* COK: +0.2000
* OCK: -0.2500

#### Pair Coeffs



### Oxirene

Structure with custom types: `HCO-COO=1(-OCO-)1=COO-HCO`
Has a `COC` cyle, with double bond between the carbons.

#### Atom Types

* COO: 87/47
* HCO: 89/46
* OCO: 122/20

#### Charges

* COO: +0.1350 \* adjusted
* HCO: +0.1150
* OCO: -0.4000


### Ethynol

Structure with custom types: `HCE-CHE#COE-OHE-HOE`

#### Atom Types

* CHE: 755/19
* COE: 759/19
* HCE: 756/46
* HOE: 97/7
* OHE: 96/5

#### Charges

* CHE: -0.2100
* COE: +0.2750 \* adjusted
* HCE: +0.2000
* HOE: +0.4180
* OHE: -0.6830
