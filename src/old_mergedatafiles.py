# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 11:35:19 2014

@author: mdoig
"""

# Program to combine data files - in order for this to work all data files
# should have existing force field parameters present. While this could lead to 
# duplicates, it is much easier and more efficient to code

# usage: python combineDatafiles.py -n 3 data1 min/max z data2 min/max z data3 min/maxz -o data.combined.out
# where -n # is the number of data files
"""
sysargv[0] : python program name (combineDatafiles.py)
sysargv[1] : -n
sysargv[2] : data1
sysargv[3] : min/max
sysargv[4] : #
sysargv[5] : data2
sysargv[6] : min/max
sysargv[7] : #
...
sysargv[] : -o
sysargv[] : outputfilename

"""
"""
Steps:
1) Read in all data files

"""

import random
import sys
import math, gzip
from commondata import readAll, getNatoms, readTS, getTSrange, getAtomType, getAtomData
from commondata import getTS, getAt, getMol, getTotMass, getCOMts, readcombmasses, readAllGzip
from commondata import getAllAtomData
from itertools import islice
import subprocess 

numdatafiles = int(sys.argv[2])
readsyntax = {}
datafile = {}
outputdata = {}
for i in range (1, numdatafiles+1):
    print i
    readsyntax[i] = {}
    readsyntax[i]['name'] = {}
    readsyntax[i]['name'] = sys.argv[(3*i)]
    readsyntax[i]['maxmin'] = sys.argv[(3*i+1 )]
    readsyntax[i]['shift'] = float(sys.argv[(3*i +2 )])
    datafile[i] = {}
    datafile[i]['atomdata'] = {}
    datafile[i]['bonddata'] = {}
    datafile[i]['angledata'] = {}
    datafile[i]['dhdata'] = {}
    datafile[i]['impdata'] = {}
    datafile[i]['masses'] = {}
    datafile[i]['boxsize'] = {}
    datafile[i]['numdata'] = {}
outputfile = sys.argv[(3*i + 4)]


#fndump = "dump.prod.all.lammpstrj"
#fndata = "datafile"

header = 9

zmax = {}
zmin = {}

for i in range(1,numdatafiles+1):
    zlist = []
    zmax[i] = {}
    zmin[i] = {}
    fndata = readsyntax[i]['name']
    print "Data file "+str(fndata)
    atomnames = getAtomType(fndata)
    print atomnames
    atomdata, bonddata, angledata, dhdata, impdata, masses, boxsize, numdata = getAllAtomData(fndata, atomnames)
    datafile[i]['atomdata'] = atomdata
    datafile[i]['bonddata'] = bonddata
    datafile[i]['angledata'] = angledata
    datafile[i]['dhdata'] = dhdata
    datafile[i]['impdata'] = impdata
    datafile[i]['masses'] = masses
    datafile[i]['boxsize'] = boxsize
    datafile[i]['numdata'] = numdata
    a = datafile[i]['atomdata'].keys()
    for j in a:
        zlist.append(float(datafile[i]['atomdata'][j]['z']))#
        
    """
    SHIFT COORDINATES BASED ON MAX/MIN
    """
    if readsyntax[i]['maxmin'] == "max":
        zmax[i] = max(zlist)
        if zmax[i] > readsyntax[i]['shift']:
            for k in a:
                atomdata[k]['z'] -= (zmax[i] - readsyntax[i]['shift'])
    if readsyntax[i]['maxmin'] == "min":
        zmin[i] = min(zlist)
        if zmin[i] < readsyntax[i]['shift']:
            for k in a:
                atomdata[k]['z'] += (readsyntax[i]['shift']- zmin[i])
    
#########################################################################################


#########################################################################################
# MERGING FILES    
# CALCULATE RENUMBERING
# #######################################################################################
outdatafile = {}
outdatafile['atomdata'] = {}
outdatafile['bonddata'] = {}
outdatafile['angledata'] = {}
outdatafile['dhdata'] = {}
outdatafile['impdata'] = {}
outdatafile['masses'] = {}
outdatafile['boxsize'] = {}
outdatafile['numdata'] = {}

atshift = 0
bshift = 0
angshift = 0
dhshift = 0
impshift = 0
attypeshift = 0
btypeshift = 0
angtypeshift = 0
dhtypeshift = 0
imptypeshift = 0

molshift = {}
moleshift = 0

for i in range(1,numdatafiles+1):
    molshift[i] = {}
    mollist = []
    l_atoms = []
    l_atoms = datafile[i]['atomdata'].keys()
    for b in l_atoms:
        mollist.append(datafile[i]['atomdata'][b]['mol'])
    molshift[i] = max(mollist)
        


for i in range (1, numdatafiles+1):
    if i == 1:
        outdatafile['atomdata'] = datafile[i]['atomdata']
    else:
        # GET LIST OF ORIGINAL KEYS
        l_atoms_orig = []
        l_atoms_orig = datafile[i]['atomdata'].keys()
        l_bonds_orig = []
        l_bonds_orig = datafile[i]['bonddata'].keys()
        l_angles_orig = []
        l_angles_orig = datafile[i]['angledata'].keys()
        l_dh_orig = []
        l_dh_orig = datafile[i]['dhdata'].keys()
        l_imp_orig = []
        l_imp_orig = datafile[i]['impdata'].keys()        
        # CALCULATE SHIFT (continuous sum over loop)
        atshift += datafile[i-1]['numdata']['natoms']
        print "atom shift - "+str(i)+" - "+str(atshift)
        attypeshift += datafile[i-1]['numdata']['natomtypes']
        print "atom type shift - "+str(i)+" - "+str(attypeshift)
        bshift += datafile[i-1]['numdata']['nbonds']
        btypeshift += datafile[i-1]['numdata']['nbondtypes']
        angshift += datafile[i-1]['numdata']['nangles']
        angtypeshift += datafile[i-1]['numdata']['nangletypes']
        dhshift += datafile[i-1]['numdata']['ndh']
        dhtypeshift += datafile[i-1]['numdata']['ndhtypes']
        impshift += datafile[i-1]['numdata']['nimp']
        imptypeshift += datafile[i-1]['numdata']['nimptypes']
        #natprev = atshift #datafile[i-1]['numdata']['natoms']
        moleshift += molshift[i-1]
        print 'moleshift '+str(moleshift)
        l_atoms = []
        l_bonds = []
        l_atoms = [x+atshift for x in l_atoms_orig]
        l_bonds = [x+bshift for x in l_bonds_orig]
        l_angles = [x+angshift for x in l_angles_orig]
        l_dh = [x+dhshift for x in l_dh_orig]
        l_imp = [x+impshift for x in l_imp_orig]
        
        #newdatafile = {}
        #newdatafile = datafile
        counter = 0
        for j, x in enumerate(l_atoms_orig):
            counter += 1
            # INCREMENT MOL TYPE + MOL NUMBER
            # TO DO
            ################################
            #datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'][l_atoms_orig[j]]

            #datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'].pop(x)
            datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'][l_atoms_orig[j]]
            #datafile[i]['atomdata'].pop(l_atoms_orig[j], None)
            datafile[i]['atomdata'][l_atoms[j]]['mol'] += moleshift
            datafile[i]['atomdata'][l_atoms[j]]['type'] += attypeshift
            print datafile[i]['atomdata'][l_atoms[j]]
            datafile[i]['atomdata'].pop(l_atoms_orig[j], None)
            
        """    
            datafile[i]['atomdata'][x]['type'] += attypeshift
            
            if j < 5:
                print "MOLSHIFT BEFORE: "+str(datafile[i]['atomdata'][x]['mol'])
            datafile[i]['atomdata'][l_atoms_orig[j]]['mol'] += moleshift
            if j < 5:
                print "MOLSHIFT AFTER: "+str(datafile[i]['atomdata'][x]['mol'])
            ####################################
            #Renumber keys in dictionary for combined atoms
            datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'][l_atoms_orig[j]]
            #if j < 5:#11750:
            print "After renumbering: "+str(datafile[i]['atomdata'][l_atoms[j]]['mol'])+" "+str(l_atoms[j])
            #print "j: "+str(j)+" x: "+str(x)+" MOLA: "+str(datafile[i]['atomdata'][l_atoms[j]]['mol'])
            #if j < 11755:
            #    if j > 11750:
            #        print "Previous mol num: "+str(datafile[i]['atomdata'][l_atoms[j-1]]['mol'])+" "+str(l_atoms[j-1])
#            print "l_atoms l_atoms_orig "+str(l_atoms[j])+" "+str(l_atoms_orig[j])
#            datafile[i]['atomdata'].pop(l_atoms_orig[j], None)              
        """
        print "Length of dictionary before pop: "+str(len(datafile[i]['atomdata']))
        """              
        for k in range (1, atshift+1):
            datafile[i]['atomdata'].pop(k, None)
        print "Length of dictionary after pop: "+str(len(datafile[i]['atomdata']))
        
        j = 0        
        """    
        """    
        for j, x in enumerate(l_atoms_orig):
             if j < 5: #11750:
                print "After renumbering part 2: "+str(datafile[i]['atomdata'][l_atoms[j]]['mol'])+" "+str(l_atoms[j])
        """
            #if j < 5:
            #    print datafile[i]['atomdata'][l_atoms[j]]['mol']
            
        #for j, x in enumerate(l_atoms_prev):
        #    newdatafile[i]['atomdata'].pop(l_atoms_orig[j], None)
            #datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'].pop(x)

        for j, x in enumerate(l_bonds_orig):
            # INCREMENT ATOM NUMBERS / BONDTYPE IN BOND LIST
            datafile[i]['bonddata'][x]['at1'] += atshift
            datafile[i]['bonddata'][x]['at2'] += atshift
            datafile[i]['bonddata'][x]['bondtype'] += btypeshift
            # SHIFT BOND NUMBERS IN DICTIONARY
            datafile[i]['bonddata'][l_bonds[j]] = datafile[i]['bonddata'].pop(x)
        #for k in range (1, bshift+1):
        #    datafile[i]['bonddata'].pop(k, None)
            
            
        for j, x in enumerate(l_angles_orig):
            # INCREMENT ATOM NUMBERS / ANGLETYPE IN ANGLE LIST
            datafile[i]['angledata'][x]['at1'] += atshift
            datafile[i]['angledata'][x]['at2'] += atshift
            datafile[i]['angledata'][x]['at3'] += atshift
            datafile[i]['angledata'][x]['angletype'] += angtypeshift
            # SHIFT ANGLE NUMBERS IN DICTIONARY
        #for k in range (1, angshift+1):
        #    datafile[i]['angledata'].pop(k, None)
            datafile[i]['angledata'][l_angles[j]] = datafile[i]['angledata'].pop(x)
            
        for j, x in enumerate(l_dh_orig):
            # INCREMENT ATOM NUMBERS / ANGLETYPE IN ANGLE LIST
            datafile[i]['dhdata'][x]['at1'] += atshift
            datafile[i]['dhdata'][x]['at2'] += atshift
            datafile[i]['dhdata'][x]['at3'] += atshift
            datafile[i]['dhdata'][x]['at4'] += atshift            
            datafile[i]['dhdata'][x]['dhtype'] += dhtypeshift
            # SHIFT ANGLE NUMBERS IN DICTIONARY
        #for k in range (1, dhshift+1):
            #datafile[i]['dhdata'].pop(k, None)
            datafile[i]['dhdata'][l_dh[j]] = datafile[i]['dhdata'].pop(x)
            
        for j, x in enumerate(l_imp_orig):
            # INCREMENT ATOM NUMBERS / ANGLETYPE IN ANGLE LIST
            datafile[i]['impdata'][x]['at1'] += atshift
            datafile[i]['impdata'][x]['at2'] += atshift
            datafile[i]['impdata'][x]['at3'] += atshift
            datafile[i]['impdata'][x]['at4'] += atshift            
            datafile[i]['impdata'][x]['imptype'] += imptypeshift
            # SHIFT ANGLE NUMBERS IN DICTIONARY
            #for k in range (1, impshift+1):
            #    datafile[i]['impdata'].pop(k, None)
            datafile[i]['impdata'][l_imp[j]] = datafile[i]['impdata'].pop(x)
            

###################################################
# WRITE OUTPUT DATA FILE

thefile = open(outputfile, "w")

###############################################
# WRITE HEADER
thefile.write("Combined LAMMPS data file\n")

totnatoms = max(datafile[numdatafiles]['atomdata'].keys())
print str(totnatoms)+" atoms"
thefile.write(str(totnatoms)+" atoms\n")
##########################################
#  NUMBERS
try:
    totnbonds = max(datafile[numdatafiles]['bonddata'].keys())
except (ValueError, TypeError):
    totnbonds = 0
print str(totnbonds)+" bonds"
thefile.write(str(totnbonds)+" bonds\n")

try:
    totnangles = max(datafile[numdatafiles]['angledata'].keys())
except (ValueError, TypeError):
    totnangles = 0
print str(totnangles)+" angles"
thefile.write(str(totnangles)+" angles\n")

try:
    totndh = max(datafile[numdatafiles]['dhdata'].keys())
except (ValueError, TypeError):
    totndh = 0
print str(totndh)+" dihedrals"
thefile.write(str(totndh)+" dihedrals\n")

try:
    totnimp = max(datafile[numdatafiles]['impdata'].keys())
except (ValueError, TypeError):
    totnimp = 0
print str(totnimp)+" impropers"
thefile.write(str(totnimp)+" impropers\n")

##############################
# TYPES
try:
    l = []
    b = datafile[numdatafiles]['atomdata'].keys()
    for i in b:
        l.append(datafile[numdatafiles]['atomdata'][i]['type'])
    totnatomstype = max(l)
except (ValueError, TypeError):
    totnatomstype = 0
print str(totnatomstype)+" atom types"
thefile.write(str(totnatomstype)+" atom types\n")

try:
    l = []
    b = datafile[numdatafiles]['bonddata'].keys()
    for i in b:
        l.append(datafile[numdatafiles]['bonddata'][i]['bondtype'])
    totnbondstype = max(l)
except (ValueError, TypeError):
    totnbondstype = 0
print str(totnbondstype)+" bond types"
thefile.write(str(totnbondstype)+" bond types\n")

try:
    l = []
    b = datafile[numdatafiles]['angledata'].keys()
    for i in b:
        l.append(datafile[numdatafiles]['angledata'][i]['angletype'])
    totnanglestype = max(l)
except (ValueError, TypeError):
    totnanglestype = 0
print str(totnanglestype)+" angle types"
thefile.write(str(totnanglestype)+" angle types\n")

try:
    l = []
    b = datafile[numdatafiles]['dhdata'].keys()
    for i in b:
        l.append(datafile[numdatafiles]['dhdata'][i]['dhtype'])
    totndhtype = max(l)
except (ValueError, TypeError):
    totndhtype = 0
print str(totndhtype)+" dihedral types"
thefile.write(str(totndhtype)+" dihedral types\n")

try:
    l = []
    b = datafile[numdatafiles]['impdata'].keys()
    for i in b:
        l.append(datafile[numdatafiles]['impdata'][i]['imptype'])
    totnimptype = max(l)
except (ValueError, TypeError):
    totnimptype = 0
print str(totnimptype)+" improper types"
thefile.write(str(totnimptype)+" improper types\n")

thefile.write("0 0 xlo xhi\n")
thefile.write("0 0 ylo yhi\n")
thefile.write("0 0 zlo zhi\n")

thefile.close()

#############################################################
#            HEADER COMPLETE
#############################################################            
            
##############################################################
#            GET X/Y/Z HI LO
            
#  TO BE COMPLETED
            
##############################################################

##############################################################
# ALL COEFFS
            
# TO BE COMPLETED
            
################################################################



###############################################################
# WRITE ATOMS PART OF DATA FILE
def writeAtomsData(outputfn, atomdata, i):
    thefile = open(outputfn,"a")
    if(i==1):
        thefile.write("\n")
        thefile.write("Atoms\n")
        thefile.write("\n")
    b = sorted(atomdata.keys())
    for i in b:
        thefile.write("%i %i %i %f %f %f %f %s %s \n" % (i, atomdata[i]['mol'], atomdata[i]['type'], \
        atomdata[i]['charge'], atomdata[i]['x'], atomdata[i]['y'], atomdata[i]['z'], '#', atomdata[i]['atomname']))
    thefile.close()
    
def writeBondData(outputfn, bonddata, i, totnbonds):
    thefile = open(outputfn,"a")
    if(i==1 and totnbonds != 0):
        thefile.write("\n")
        thefile.write("Bonds\n")
        thefile.write("\n")
    if(totnbonds != 0):
        b = sorted(bonddata.keys())
        for i in b:
            thefile.write("%i %i %i %i \n" % (i, bonddata[i]['bondtype'], \
            bonddata[i]['at1'], bonddata[i]['at2']))
    thefile.close()

def writeAngleData(outputfn, angledata, i, totnangles):
    thefile = open(outputfn,"a")
    if(i==1 and totnangles != 0):
        thefile.write("\n")
        thefile.write("Angles\n")
        thefile.write("\n")
    if(totnangles != 0):       
        b = sorted(angledata.keys())
        for i in b:
            thefile.write("%i %i %i %i %i \n" % (i, angledata[i]['angletype'], \
            angledata[i]['at1'], angledata[i]['at2'], angledata[i]['at3']))
    thefile.close()

def writeDHData(outputfn, dhdata, i, totndh):
    thefile = open(outputfn,"a")
    if(i==1 and totndh != 0):
        thefile.write("\n")
        thefile.write("Dihedrals\n")
        thefile.write("\n")
    if(totndh != 0):       
        b = sorted(dhdata.keys())
        for i in b:
            thefile.write("%i %i %i %i %i %i \n" % (i, dhdata[i]['dhtype'], \
            dhdata[i]['at1'], dhdata[i]['at2'], dhdata[i]['at3'], dhdata[i]['at4']))
    thefile.close()

def writeIMPData(outputfn, impdata, i, totnimp):
    thefile = open(outputfn,"a")
    if(i==1 and totnimp != 0):
        thefile.write("\n")
        thefile.write("Impropers\n")
        thefile.write("\n")
    if(totnimp != 0):       
        b = sorted(impdata.keys())
        for i in b:
            thefile.write("%i %i %i %i %i %i \n" % (i, impdata[i]['imptype'], \
            impdata[i]['at1'], impdata[i]['at2'], impdata[i]['at3'], impdata[i]['at4']))
    thefile.close()

def writeMasses(outputfn, massdata, i, totnatomstype):
    thefile = open(outputfn, "a")
    if (i==1):
        thefile.write("\n")
        thefile.write("Masses\n")
        thefile.write("\n")
    b = sorted(massdata.keys())
    counter = 0
    for i, x in enumerate(b):
        counter += 1
        thefile.write("%i %f %s %s \n" % (counter, massdata[x], "#", x))


for i in range(1, numdatafiles+1):
    writeMasses(outputfile,datafile[i]['masses'], i, totnatomstype)

for i in range (1, numdatafiles+1):
    writeAtomsData(outputfile, datafile[i]['atomdata'], i)
    
for i in range(1, numdatafiles+1):
    writeBondData(outputfile, datafile[i]['bonddata'], i, totnbonds)
    
for i in range(1, numdatafiles+1):
    writeAngleData(outputfile, datafile[i]['angledata'], i, totnangles)
 
for i in range(1, numdatafiles+1):
    writeDHData(outputfile, datafile[i]['dhdata'], i, totndh)
    
for i in range(1, numdatafiles+1):
    writeIMPData(outputfile, datafile[i]['impdata'], i, totnimp)   
 