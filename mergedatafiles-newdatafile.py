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

try:
    temp = sys.argv[1]
except:
    print "Program usage: python mergedatafiles-newdatafile.py  -n 3 data1 min/max z data2 min/max z data3 min/maxz -o data.combined.out (where -n # is the number of data files)"
    sys.exit(-1)


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
    datafile[i]['paircoeff'] = {}
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
    atomdata, bonddata, angledata, dhdata, impdata, masses, boxsize, numdata, paircoeff, bondcoeff, anglecoeff, dhcoeff, impcoeff = getAllAtomData(fndata, atomnames)
    datafile[i]['atomdata'] = atomdata
    datafile[i]['bonddata'] = bonddata
    datafile[i]['angledata'] = angledata
    datafile[i]['dhdata'] = dhdata
    datafile[i]['impdata'] = impdata
    datafile[i]['masses'] = masses
    datafile[i]['boxsize'] = boxsize
    datafile[i]['numdata'] = numdata
    datafile[i]['paircoeff'] = paircoeff
    datafile[i]['bondcoeff'] = bondcoeff
    datafile[i]['anglecoeff'] = anglecoeff
    datafile[i]['dhcoeff'] = dhcoeff
    datafile[i]['impcoeff'] = impcoeff
    
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
    #print "MOLLIST - EMPTY:"
    #print mollist
    l_atoms = []
    l_atoms = datafile[i]['atomdata'].keys()
    for b in l_atoms:
        mollist.append(datafile[i]['atomdata'][b]['mol'])
    molshift[i] = max(mollist)
    #print "MAX: MOL LIST"+str(molshift[i])
        


for i in range (1, numdatafiles+1):
    if i == 1:
        outdatafile['atomdata'] = datafile[i]['atomdata']
        newdatafile = {}
        newdatafile[i] = {}
        newdatafile[i] = datafile[i]
        
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
        
        
        newdatafile[i] = {}
        newdatafile[i]['atomdata']= {}
        newdatafile[i]['masses'] = {}
        newdatafile[i]['masses'] = datafile[i]['masses']
        newdatafile[i]['numdata'] = {}
        newdatafile[i]['numdata'] = datafile[i]['numdata']        
        counter = 0
        for j, x in enumerate(l_atoms_orig):
            counter += 1
            # INCREMENT MOL TYPE + MOL NUMBER
            # TO DO
            ################################
            #datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'][l_atoms_orig[j]]

            #datafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'].pop(x)
            newdatafile[i]['atomdata'][l_atoms[j]] = {}
            newdatafile[i]['atomdata'][l_atoms[j]] = datafile[i]['atomdata'][l_atoms_orig[j]]
            #print "OLD/NEW"
            #print datafile[i]['atomdata'][l_atoms_orig[j]]
            #print newdatafile[i]['atomdata'][l_atoms[j]]
            #datafile[i]['atomdata'].pop(l_atoms_orig[j], None)
            #print "BEFORE"
            #print newdatafile[i]['atomdata'][l_atoms[j]]
            newdatafile[i]['atomdata'][l_atoms[j]]['mol'] += moleshift
            newdatafile[i]['atomdata'][l_atoms[j]]['type'] += attypeshift
            #print "AFTER:" 
            #print newdatafile[i]['atomdata'][l_atoms[j]]
            #newdatafile[i]['atomdata'].pop(l_atoms_orig[j], None)
        
        #for j, x in enumerate(l_atoms_orig):          
        #    print "MOL TYPE - NEW: "+str(newdatafile[i]['atomdata'][l_atoms[j]]['mol'])
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
        #print "Length of dictionary before pop: "+str(len(newdatafile[i]['atomdata']))
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
        
        newdatafile[i]['bonddata'] = {}
        for j, x in enumerate(l_bonds_orig):
            # INCREMENT ATOM NUMBERS / BONDTYPE IN BOND LIST
            newdatafile[i]['bonddata'][l_bonds[j]] = {}
            newdatafile[i]['bonddata'][l_bonds[j]] = datafile[i]['bonddata'][l_bonds_orig[j]]
            newdatafile[i]['bonddata'][l_bonds[j]]['at1'] += atshift
            newdatafile[i]['bonddata'][l_bonds[j]]['at2'] += atshift
            newdatafile[i]['bonddata'][l_bonds[j]]['bondtype'] += btypeshift
            # SHIFT BOND NUMBERS IN DICTIONARY
            #newdatafile[i]['bonddata'].pop(l_bonds_orig[j], None)
            #newdatafile[i]['bonddata'][l_bonds[j]] = datafile[i]['bonddata'].pop(x)
        #for k in range (1, bshift+1):
        #    datafile[i]['bonddata'].pop(k, None)
            
        newdatafile[i]['angledata'] = {}    
        for j, x in enumerate(l_angles_orig):
            # INCREMENT ATOM NUMBERS / ANGLETYPE IN ANGLE LIST
            newdatafile[i]['angledata'][l_angles[j]] = {}
            newdatafile[i]['angledata'][l_angles[j]] = datafile[i]['angledata'][l_angles_orig[j]]
            newdatafile[i]['angledata'][l_angles[j]]['at1'] += atshift
            newdatafile[i]['angledata'][l_angles[j]]['at2'] += atshift
            newdatafile[i]['angledata'][l_angles[j]]['at3'] += atshift
            newdatafile[i]['angledata'][l_angles[j]]['angletype'] += angtypeshift           
            # SHIFT ANGLE NUMBERS IN DICTIONARY
            #newdatafile[i]['angledata'].pop(l_angles_orig[j], None)
            #datafile[i]['angledata'][l_angles[j]] = datafile[i]['angledata'].pop(x)
        
        newdatafile[i]['dhdata'] = {}
        for j, x in enumerate(l_dh_orig):
            # INCREMENT ATOM NUMBERS / ANGLETYPE IN ANGLE LIST
            newdatafile[i]['dhdata'][l_dh[j]] = {}
            newdatafile[i]['dhdata'][l_dh[j]] = datafile[i]['dhdata'][l_dh_orig[j]]
            newdatafile[i]['dhdata'][l_dh[j]]['at1'] += atshift
            newdatafile[i]['dhdata'][l_dh[j]]['at2'] += atshift
            newdatafile[i]['dhdata'][l_dh[j]]['at3'] += atshift
            newdatafile[i]['dhdata'][l_dh[j]]['at4'] += atshift
            newdatafile[i]['dhdata'][l_dh[j]]['dhtype'] += dhtypeshift      
            # SHIFT ANGLE NUMBERS IN DICTIONARY
        #for k in range (1, dhshift+1):
            #datafile[i]['dhdata'].pop(k, None)
            #datafile[i]['dhdata'][l_dh[j]] = datafile[i]['dhdata'].pop(x)
            #newdatafile[i]['dhdata'].pop(l_dh_orig[j], None)
        
        newdatafile[i]['impdata'] = {}
        for j, x in enumerate(l_imp_orig):
            # INCREMENT ATOM NUMBERS / ANGLETYPE IN ANGLE LIST
            newdatafile[i]['impdata'][l_imp[j]] = {}
            newdatafile[i]['impdata'][l_imp[j]] = datafile[i]['impdata'][l_imp_orig[j]]
            newdatafile[i]['impdata'][l_imp[j]]['at1'] += atshift
            newdatafile[i]['impdata'][l_imp[j]]['at2'] += atshift
            newdatafile[i]['impdata'][l_imp[j]]['at3'] += atshift
            newdatafile[i]['impdata'][l_imp[j]]['at4'] += atshift
            newdatafile[i]['impdata'][l_imp[j]]['imptype'] += imptypeshift   
                        # SHIFT ANGLE NUMBERS IN DICTIONARY
            #for k in range (1, impshift+1):
            #    datafile[i]['impdata'].pop(k, None)
            #newdatafile[i]['impdata'].pop(l_imp_orig[j], None)
            #datafile[i]['impdata'][l_imp[j]] = datafile[i]['impdata'].pop(x)
            

###################################################
# WRITE OUTPUT DATA FILE

thefile = open(outputfile, "w")

###############################################
# WRITE HEADER
thefile.write("Combined LAMMPS data file\n")

totnatoms = max(newdatafile[numdatafiles]['atomdata'].keys())
print str(totnatoms)+" atoms"
thefile.write(str(totnatoms)+" atoms\n")
##########################################
totnbonds = 0
totnangles = 0
totndh = 0
totnimp = 0
for i in range(1,numdatafiles+1):
    try:
        totnbonds += len(newdatafile[i]['bonddata'].keys())
        totnangles += len(newdatafile[i]['angledata'].keys())
        totndh += len(newdatafile[i]['dhdata'].keys())
        totnimp += len(newdatafile[i]['imp'].keys())
    except:
        junk = 0

#print totnbonds, totnangles, totndh, totnimp
thefile.write(str(totnbonds)+" bonds\n")
thefile.write(str(totnangles)+" angles\n")
thefile.write(str(totndh)+" dihedrals\n")
thefile.write(str(totnimp)+" impropers\n")

"""
DEPRECATED - 26/11/14
try:
    totnbonds = max(newdatafile[numdatafiles]['bonddata'].keys())
except (ValueError, TypeError):
    totnbonds = 0


try:
    totnangles = max(newdatafile[numdatafiles]['angledata'].keys())
except (ValueError, TypeError):
    totnangles = 0
print str(totnangles)+" angles"
thefile.write(str(totnangles)+" angles\n")

try:
    totndh = max(newdatafile[numdatafiles]['dhdata'].keys())
except (ValueError, TypeError):
    totndh = 0
print str(totndh)+" dihedrals"
thefile.write(str(totndh)+" dihedrals\n")

try:
    totnimp = max(newdatafile[numdatafiles]['impdata'].keys())
except (ValueError, TypeError):
    totnimp = 0
print str(totnimp)+" impropers"
thefile.write(str(totnimp)+" impropers\n")
"""

##############################
# TYPES
lat = []
lbo = []
lan = []
ldh = []
lim = []

for i in range(1,numdatafiles+1):
    try:
        at = newdatafile[i]['atomdata'].keys()
        bo = newdatafile[i]['bonddata'].keys()
        an = newdatafile[i]['angledata'].keys()
        dh = newdatafile[i]['dhdata'].keys()
        im = newdatafile[i]['impdata'].keys()
        for j in at:
            lat.append(newdatafile[i]['atomdata'][j]['type'])
        for j in bo:
            lbo.append(newdatafile[i]['bonddata'][j]['bondtype'])
        for j in an:
            lan.append(newdatafile[i]['angledata'][j]['angletype'])
        for j in dh:
            ldh.append(newdatafile[i]['dhdata'][j]['dhtype'])
        for j in im:
            lim.append(newdatafile[i]['impdata'][j]['imptype'])
    except:
        junk = 0
try:
    totnatomstype = max(lat)
except:
    totnatomstype = 0
try:
    totnbondstype = max(lbo)
except:
    totnbondstype = 0
try:
    totnanglestype = max(lan)
except:
    totnanglestype = 0
try:
    totndhtype = max(ldh)
except:
    totndhtype = 0
try:
    totnimptype = max(lim)
except:
    totnimptype = 0
#print totnbondstype, totnanglestype, totndhtype, totnimptype
           
print str(totnatomstype)+" atom types"           
print str(totnbondstype)+" bond types"
print str(totnanglestype)+" angle types"
print str(totndhtype)+" dihedral types"
print str(totnimptype)+" improper types"           

thefile.write(str(totnatomstype)+" atom types\n")
thefile.write(str(totnbondstype)+" bond types\n")
thefile.write(str(totnanglestype)+" angle types\n")
thefile.write(str(totndhtype)+" dihedral types\n")
thefile.write(str(totnimptype)+" improper types\n")
thefile.write("0 0 xlo xhi\n")
thefile.write("0 0 ylo yhi\n")
thefile.write("0 0 zlo zhi\n")

thefile.close()


"""
DEPRECATED - 26/11/14
try:
    l = []
    b = newdatafile[numdatafiles]['atomdata'].keys()
    for i in b:
        l.append(newdatafile[numdatafiles]['atomdata'][i]['type'])
    totnatomstype = max(l)
except (ValueError, TypeError):
    totnatomstype = 0
print str(totnatomstype)+" atom types"


try:
    l = []
    b = newdatafile[numdatafiles]['bonddata'].keys()
    for i in b:
        l.append(newdatafile[numdatafiles]['bonddata'][i]['bondtype'])
    totnbondstype = max(l)
except (ValueError, TypeError):
    totnbondstype = 0
print str(totnbondstype)+" bond types"
thefile.write(str(totnbondstype)+" bond types\n")

try:
    l = []
    b = newdatafile[numdatafiles]['angledata'].keys()
    for i in b:
        l.append(newdatafile[numdatafiles]['angledata'][i]['angletype'])
    totnanglestype = max(l)
except (ValueError, TypeError):
    totnanglestype = 0
print str(totnanglestype)+" angle types"
thefile.write(str(totnanglestype)+" angle types\n")

try:
    l = []
    b = newdatafile[numdatafiles]['dhdata'].keys()
    for i in b:
        l.append(newdatafile[numdatafiles]['dhdata'][i]['dhtype'])
    totndhtype = max(l)
except (ValueError, TypeError):
    totndhtype = 0
print str(totndhtype)+" dihedral types"
thefile.write(str(totndhtype)+" dihedral types\n")

try:
    l = []
    b = newdatafile[numdatafiles]['impdata'].keys()
    for i in b:
        l.append(newdatafile[numdatafiles]['impdata'][i]['imptype'])
    totnimptype = max(l)
except (ValueError, TypeError):
    totnimptype = 0
print str(totnimptype)+" improper types"
thefile.write(str(totnimptype)+" improper types\n")
"""


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

def writeMasses(outputfn, massdata, i, totnatomstype, counter):
    thefile = open(outputfn, "a")
    if (i==1):
        thefile.write("\n")
        thefile.write("Masses\n")
        thefile.write("\n")
    b = sorted(massdata.keys())
    for i, x in enumerate(b):
        thefile.write("%i %f %s %s \n" % (counter, massdata[x], "#", x))
        counter += 1
    return counter

def writePairCoeffs(outputfn, paircoeffs, i, totnatomstype, counter):
    thefile = open(outputfn, "a")
    if (i==1):
        thefile.write("\n")
        thefile.write("Pair Coeffs\n")
        thefile.write("\n")
    b = sorted(paircoeffs.keys())
    #print "B - sorted pair coeffs:"
    #print b
    #print "pair coeffs"
    #print paircoeffs
    #print "pair coeff"
    #print paircoeff
    for i in b:
        #print "i (in b) :"+str(i)
        #print "counter :"+str(counter)
        #print "SAME AS FILE WRITE"
        #print counter, paircoeff[i]['eps'], paircoeff[i]['sig'], '#', paircoeff[i]['label']
        thefile.write("%i %s %s %s %s \n" % (counter, paircoeffs[i]['eps'], paircoeffs[i]['sig'], '#', paircoeffs[i]['label']))
        counter += 1
    thefile.close()
    return counter
    
def writeBondCoeffs(outputfn, bondcoeffs, i, totnbondstype, counter):
    thefile = open(outputfn, "a")
    if (i==1):
        thefile.write("\n")
        thefile.write("Bond Coeffs\n")
        thefile.write("\n")
    b = sorted(bondcoeffs.keys())
    for i in b:
        thefile.write("%i %s %s %s %s \n" % (counter, bondcoeffs[i]['k'], bondcoeffs[i]['r'], '#', bondcoeffs[i]['label']))
        counter += 1
    thefile.close()
    return counter
    
def writeAngleCoeffs(outputfn, anglecoeffs, i, totnanglestype, counter):
    thefile = open(outputfn, "a")
    if (i==1):
        thefile.write("\n")
        thefile.write("Angle Coeffs\n")
        thefile.write("\n")
    b = sorted(anglecoeffs.keys())
    for i in b:
        thefile.write("%i %s %s %s %s \n" % (counter, anglecoeffs[i]['k'], anglecoeffs[i]['theta'], '#', anglecoeffs[i]['label']))
        counter += 1
    thefile.close()
    return counter
    
def writeDHCoeffs(outputfn, dhcoeffs, i, totndhtype, counter):
    thefile = open(outputfn, "a")
    if (i==1):
        thefile.write("\n")
        thefile.write("Dihedral Coeffs\n")
        thefile.write("\n")
    b = sorted(dhcoeffs.keys())
    for i in b:
        thefile.write("%i %s %s %s %s %s %s \n" % (counter, dhcoeffs[i]['a'], dhcoeffs[i]['b'], dhcoeffs[i]['c'], dhcoeffs[i]['d'], '#', dhcoeffs[i]['label']))
        counter += 1
    thefile.close()
    return counter

def writeIMPCoeffs(outputfn, impcoeffs, i, totnimptype, counter):
    thefile = open(outputfn, "a")
    if (i==1 and totnimp != 0):
        thefile.write("\n")
        thefile.write("Improper Coeffs\n")
        thefile.write("\n")
    b = sorted(impcoeffs.keys())
    for i in b:
        thefile.write("%i %s %s %s %s %s %s \n" % (counter, impcoeffs[i]['a'], impcoeffs[i]['b'], impcoeffs[i]['c'], impcoeffs[i]['d'], '#', impcoeffs[i]['label']))
        counter += 1
    thefile.close()
    return counter

paircount = 1
for i in range(1, numdatafiles+1):
    paircount = writePairCoeffs(outputfile, datafile[i]['paircoeff'], i , totnatomstype, paircount)

bondcount = 1
for i in range(1, numdatafiles+1):
    bondcount = writeBondCoeffs(outputfile, datafile[i]['bondcoeff'], i , totnbondstype, bondcount)

anglecount = 1
for i in range(1, numdatafiles+1):
    anglecount = writeAngleCoeffs(outputfile, datafile[i]['anglecoeff'], i , totnanglestype, anglecount)

dhcount = 1
for i in range(1, numdatafiles+1):
    dhcount = writeDHCoeffs(outputfile, datafile[i]['dhcoeff'], i , totndhtype, dhcount)
    
impcount = 1
for i in range(1, numdatafiles+1):
    impcount = writeIMPCoeffs(outputfile, datafile[i]['impcoeff'], i , totnimptype, impcount)

masscount = 1
for i in range(1, numdatafiles+1):
    masscount = writeMasses(outputfile,newdatafile[i]['masses'], i, totnatomstype, masscount)
    
for i in range (1, numdatafiles+1):
    writeAtomsData(outputfile, newdatafile[i]['atomdata'], i)
    
for i in range(1, numdatafiles+1):
    writeBondData(outputfile, newdatafile[i]['bonddata'], i, totnbonds)
    
for i in range(1, numdatafiles+1):
    writeAngleData(outputfile, newdatafile[i]['angledata'], i, totnangles)
 
for i in range(1, numdatafiles+1):
    writeDHData(outputfile, newdatafile[i]['dhdata'], i, totndh)
    
for i in range(1, numdatafiles+1):
    writeIMPData(outputfile, newdatafile[i]['impdata'], i, totnimp)   
 