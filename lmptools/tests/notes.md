Test molecule generation and configuration
==========================================

Three molecules were generated, for testing purposes: ketene (1), oxirene (2), and ethynol (3).
They were chosen because they have the same number of atoms, but offer distinct geometries and number of atom types.
Values are not guaranteed to be accurate, but are sufficient for testing purposes.
Steps to generate all three molecules are shown below, with commands substituted with `<variables>` where appliable, and `<variable>` values are detailed elsewhere.

NOTE: should these molecules be recreated, the test values will need to be updated to match the newly generated coordinates.

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
which contains the charge values for each atom, then set the pair, bond, angle, dihedral, and improper parameters, see the corresponding section for each molecule, below.

`.lammps` files were then compressed to `.gz` using the command:
```
gzip -k <molecule>.lammps
```

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

#### OPLS parameters

```
 Pair Coeffs

1 0.086 3.3  # 1  CHK 900/47
2 0.086 3.3  # 2  COK 904/110
3 0.03  2.42 # 3  HCK 899/46
4 0.21  2.96 # 4  OCK 905/4

 Bond Coeffs

1 350 1.305 # 1  CHK-COK 47-110
2 170 1.08  # 2  CHK-HCK 46-47
3 350 1.171 # 3  COK-OCK 4-110

 Angle Coeffs

1 80.0 180.0 # 1  CHK-COK-OCK 4-110-47
2 17.5 117.0 # 2  COK-CHK-HCK 46-47-46
4 17.5 117.0 # 3  HCK-CHK-HCK 46-47-46

 Dihedral Coeffs

1 0.0 0.0 0.0 0.0 # 1  HCK-CHK-COK-OCK 46-47-110-4

 Improper Coeffs

1 15.0 180.0 # 1  COK-HCK-CHK-HCK 110-46-47-46
```

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

#### OPLS parameters

```
 Pair Coeffs

1 0.076 3.75 # 1  COO 87/47
2 0.03  2.42 # 2  HCO 89/46
3 0.14  2.9  # 3  OCO 122/20

 Bond Coeffs

1 274.5 1.34 # 1  COO-COO 47-47
2 170.0 1.08 # 2  COO-HCO 46-47
3 225.0 1.37 # 3  COO-OCO 20-47

 Angle Coeffs

1 17.5 120.0 # 1  COO-COO-HCO 46-47-47
2 35.0 123.0 # 2  COO-COO-OCO 20-47-47
3 37.5  60.0 # 3  COO-OCO-COO 47-20-47 modified to smaller angle
4 17.5 114.5 # 4  HCO-COO-OCO 20-47-46

 Dihedral Coeffs

1 -3.5  3.0 0.0  0.0 # 1  COO-COO-OCO-COO 47-20-47-47 13-20-47-47
2  0.0 14.0 0.0  0.0 # 2  HCO-COO-COO-HCO 46-47-47-46
3  0.0 14.0 0.0  0.0 # 3  HCO-COO-COO-OCO 20-47-47-46
4  0.0  0.0 0.76 0.0 # 4  HCO-COO-OCO-COO 47-20-47-46 13-20-47-36
5  0.0 14.0 0.0  0.0 # 5  OCO-COO-COO-OCO 20-47-47-20 0-47-47-0

 Improper Coeffs

1 15.0 180.0 # 1  COO-HCO-COO-OCO 47-46-47-20
```

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

#### OPLS parameters

```
 Pair Coeffs

1 0.086 3.3  # 1  CHE 755/19
2 0.1   3.3  # 2  COE 759/19
3 0.015 2.42 # 3  HCE 756/46
4 0.0   0.0  # 4  HOE 97/7
5 0.17  3.12 # 5  OHE 96/5

 Bond Coeffs

1 575.0 1.21  # 1  CHE-COE 19-19
2 210.0 1.08  # 2  CHE-HCE 19-46
3 275.0 0.137 # 3  COE-OHE 5-19 5-47
4 276.5 0.945 # 4  HOE-OHE 5-7

 Angle Coeffs

1 75.0 180.0 # 1  CHE-COE-OHE 5-19-19 13-19-19
2 56.0 180.0 # 2  COE-CHE-HCE 19-19-46
3 17.5 109.0 # 3  COE-OHE-HOE 7-5-19 7-5-47

 Dihedral Coeffs

1 0.0 0.0 0.0 0.0 # 1  CHE-COE-OHE-HOE 7-5-19-19 0-19-19-0
2 0.0 0.0 0.0 0.0 # 2  HCE-CHE-COE-OHE 5-19-19-46 0-19-19-0
```
