########

#module containing commonly used data reading and analysis tools

# import command - to be pasted in program code
"""
from commondata import readAll, getNatoms, readTS, getTSrange, getAtomType, getAtomData
from commondata import getTS, getAt, getMol, getTotMass, getCOMts, readcombmasses, readAllGzip
"""
###################################################################
###################################################################
"""
List of modules
1) readAll
2) getNatoms
3) readTS
4) getTSrange
5) getAtomType
6) getAtomData
7) getTS
8) getAt
9) getMol
10) getTotMass
11) getCOMTS
12) readcombmasses

"""
import gzip

#############################################################
"""
 readAll - Reads the entire file in line-by-line
 takes a filename, returns a list of line-by-line read in
 last modified: 13/03/14
 """
def readAll(fn):
    file = open(fn, "r")
    lines = file.readlines()
    file.close()
    return lines
#############################################################


"""
 readAllGzip - Reads the entire file in line-by-line
 takes a filename, returns a list of line-by-line read in
 last modified: 22/04/14
 """
def readAllGzip(fn):
    file = gzip.open(fn, "rb")
    lines = file.readlines()
    file.close()
    return lines
#############################################################

"""
getNatoms - gets the total number of atoms in dump file
last modified: 13/03/14
"""
def getNatoms(lines):
    for i, x in enumerate(lines):
        if lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i+1])
            return natoms

############################################################

"""
getTSrange
last modified: 13/03/14
"""
def getTSrange(lines):
    tsrange = {}
    tscount = 1
    for i,x in enumerate(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1].strip("\n"))
            tsrange[tscount] = ts
            tscount += 1
    return tsrange


"""
getTS
last modified: 13/03/14
"""
# select timestep range
# ts1 = first time step or 'first''
# ts2 = second time step or 'last'
def getTS(traj, atomdata, ts1, ts2):
    if ts1 == 'first':
        ts1 = sorted(traj.keys())[0]
    if ts2 == 'last':
        ts2 = sorted(traj.keys())[len(traj)-1]
    sortedts = sorted(traj.keys())
    tsrange = []
    for i, x in enumerate(sortedts):
        if int(x) >= int(ts1):
            if int(x) <= int(ts2):
                tsrange.append(x)
    return tsrange

###########################################################
"""
readTS
last modified: 13/03/14
"""
def readTS(lines, header, natoms, tsnum, atomnames):
    traj = {}
    j = 0
    k = 0
    nlines = natoms + header
    firstline = ((tsnum-1)*nlines)
    lastline = firstline + nlines-1
    #print 'firstline/lastline '+str(firstline)+' / '+str(lastline)
    #for i, x in enumerate(lines[firstline:lastline]):
    for i in range (firstline, lastline):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1].strip("\n"))
            traj[ts] = {}
            traj[ts]["atom"] = {}
            traj[ts]["boxsize"] = {}
            """
            EDITS 23/06/15
            """
            traj[ts]["boxx"] = {}
            traj[ts]["boxy"] = {}
            traj[ts]["boxz"] = {}
            """
            END OF EDITS 23/06/15
            """
        if "BOX BOUNDS" in lines[i]:
            boxx = map(float,lines[i+1].split())
            boxy = map(float, lines[i+2].split())
            boxz = map(float, lines[i+3].split())
            header = 9
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            """
            EDITS 23/06/15
            """
            traj[ts]["boxx"] = (boxx)
            traj[ts]["boxy"] = (boxy)
            traj[ts]["boxz"] = (boxz)
            """
            END OF EDITS 23/06/15
            """
            #traj[ts]["boxsize"]['lx'] = lx
            #traj[ts]["boxsize"]['ly'] = ly
            #traj[ts]["boxsize"]['lz'] = lz
            traj[ts]["boxsize"] = (lx, ly, lz)
        if "xy xz yz pp pp" in lines[i]:
            boxx = map(float,lines[i+1].split())
            boxy = map(float, lines[i+2].split())
            boxz = map(float, lines[i+3].split())
            header = 9
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            """
            EDITS 23/06/15
            """
            traj[ts]["boxx"] = (boxx)
            traj[ts]["boxy"] = (boxy)
            traj[ts]["boxz"] = (boxz)
            """
            END OF EDITS 23/06/15
            """
            traj[ts]["boxsize"] = (lx, ly, lz)
        if lines[i].startswith("ITEM: ATOMS"):
            elems = lines[i].split()
            for j in range(natoms):
                el = lines[i+j+1].split()
                for l, k in enumerate(elems[2:]):
                    if l==0:
                        atomid=int(el[l])
                        traj[ts]["atom"][atomid]={}
                        continue
                    traj[ts]["atom"][atomid][k] = el[l]
                #Change from numerical atom type (from dump) to atom type name (from data file)
                traj[ts]["atom"][atomid]["type"] = atomnames[traj[ts]["atom"][atomid]["type"]]
    return traj







#########################################################
"""
getAtomType
reads data file (fn) and return a list of atom numbers with all spaces/new line characters removed
 ln = (line.split("#")[1]) splits at the '#' character and puts everything after # into ln
 .replace(" ","") replaces all spaces with no spaces
 .strip("\n") removes any new line characters
 last modified: 13/03/14
"""
def getAtomType(fn):
    file = open(fn, "r")
    line = ""
    junk = ""
    while True:
        line = file.readline()
        if "atom types" in line:
            sp = line.split()
            print "split at atom types "
            natomtypes = int(sp[0])
            #print 'natom types: '+str(natomtypes)
            atomnames = {}
        if "Pair Coeffs" in line:
            junk = file.readline() # Reads blank line
            for i in range(natomtypes):
                line = file.readline()
                #print "line: "+str(line)
                ln = ((line.split("#")[1]).replace(" ","")).strip("\n")
                atomnames[str(i+1)] = ln
                #atomnames.append(ln.strip("\n"))
            file.close()
            return atomnames
#############################################################



# Not done from here down:


"""
getAtomData
last modified: 13/03/14
"""
def getAtomData(fn, atomnames):
    lines = readAll(fn)
    junk = ""
    atomdata = {}
    mass = {}
    for i, x in enumerate(lines):
        if "atoms" in lines[i]:
            natoms = ((lines[i].split("atoms")[0]).replace(" ","")).strip("\n")

        if "atom types" in lines[i]:
            ntypes = ((lines[i].split("atom")[0]).replace(" ","")).strip("\n")
            print "ntypes "+str(ntypes)

        if "Masses" in lines[i]:
            junk = lines[i+1]
            masscol = ['mass', 'type']
            for j in range(int(ntypes)):
                el = lines[i+2+j].split()
                print el
                mass[el[3]] = float(el[1])
                #print mass[el[4]]

        if "Atoms" in lines[i]:
            junk = lines[i+1]
            columns = ['mol', 'type', 'charge', 'x', 'y', 'z']
            print columns
            print columns[1]
            for j in range(int(natoms)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        atomid = int(el[0])
                        atomdata[atomid]= {}
                        atomdata[atomid][str(k)] = {}
                    atomdata[atomid][str(k)] = el[l+1]
                    #print atomnames
                atomdata[atomid]['type'] = atomnames[atomdata[atomid]['type']]
    return atomdata, mass

##########################################################

def getAllAtomData(fn, atomnames):
    lines = readAll(fn)
    junk = ""
    atomdata = {}
    bonddata = {}
    angledata = {}
    dhdata = {}
    impdata = {}
    mass = {}
    boxsize = {}
    numdata = {}
    pair = {}
    bond = {}
    angle = {}
    dh = {}
    imp = {}

    numdata['natoms'] = {}
    numdata['nbonds'] = {}
    numdata['nangles'] = {}
    numdata['ndh'] = {}
    numdata['nimp'] = {}

    numdata['natomtypes'] = {}
    numdata['nbondtypes'] = {}
    numdata['nangletypes'] = {}
    numdata['ndhtypes'] = {}
    numdata['nimptypes'] = {}

    for i, x in enumerate(lines):
        if "atoms" in lines[i]:
            natoms = int(((lines[i].split("atoms")[0]).replace(" ","")).strip("\n"))
            numdata['natoms'] = natoms
        if "bonds" in lines[i]:
            nbonds = int(((lines[i].split("bonds")[0]).replace(" ","")).strip("\n"))
            numdata['nbonds'] = nbonds
        if "angles" in lines[i]:
            nangles = int(((lines[i].split("angles")[0]).replace(" ","")).strip("\n"))
            numdata['nangles'] = nangles
        if "dihedrals" in lines[i]:
            ndh = int(((lines[i].split("dihedrals")[0]).replace(" ","")).strip("\n"))
            numdata['ndh'] = ndh
        if "impropers" in lines[i]:
            nimp = int(((lines[i].split("impropers")[0]).replace(" ","")).strip("\n"))
            numdata['nimp'] = nimp
        if "atom types" in lines[i]:
            natomtypes = int(((lines[i].split("atom")[0]).replace(" ","")).strip("\n"))
            numdata['natomtypes'] = natomtypes
            print "natomtypes "+str(natomtypes)
        if "bond types" in lines[i]:
            nbondtypes = int(((lines[i].split("bond")[0]).replace(" ","")).strip("\n"))
            numdata['nbondtypes'] = nbondtypes
            print "nbondtypes "+str(nbondtypes)
        if "angle types" in lines[i]:
            nangletypes = int(((lines[i].split("angle")[0]).replace(" ","")).strip("\n"))
            numdata['nangletypes'] = nangletypes
            print "nangletypes "+str(nangletypes)
        if "dihedral types" in lines[i]:
            ndhtypes = int(((lines[i].split("dihedral")[0]).replace(" ","")).strip("\n"))
            numdata['ndhtypes'] = ndhtypes
            print "ndhtypes "+str(ndhtypes)
        if "improper types" in lines[i]:
            nimptypes = int(((lines[i].split("improper")[0]).replace(" ","")).strip("\n"))
            numdata['nimptypes'] = nimptypes
            print "nimptypes "+str(nimptypes)
        if "xlo xhi" in lines[i]:
            ln = lines[i].split()
            xlo = float(ln[0])
            xhi = float(ln[1])
            boxsize['xlo'] = {}
            boxsize['xhi'] = {}
            boxsize['xlo'] = xlo
            boxsize['xhi'] = xhi
        if "ylo yhi" in lines[i]:
            ln = lines[i].split()
            ylo = float(ln[0])
            yhi = float(ln[1])
            boxsize['ylo'] = {}
            boxsize['yhi'] = {}
            boxsize['ylo'] = ylo
            boxsize['yhi'] = yhi
        if "zlo zhi" in lines[i]:
            ln = lines[i].split()
            zlo = float(ln[0])
            zhi = float(ln[1])
            boxsize['zlo'] = {}
            boxsize['zhi'] = {}
            boxsize['zlo'] = zlo
            boxsize['zhi'] = zhi


        if "Pair Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ 'eps', 'sig', '#', 'label']
            for j in range(int(natomtypes)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l==0:
                        paircoeffid = int(el[0])
                        pair[paircoeffid] = {}
                        pair[paircoeffid][str(k)] = {}
                    pair[paircoeffid][str(k)] = (el[l+1])
                    #pair[paircoeffid]['eps'] = float(pair[paircoeffid]['eps'])
                    #pair[paircoeffid]['sig'] = float(pair[paircoeffid]['sig'])
            print "Pair Coeffs:"
            print pair

        if "Bond Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ 'k', 'r', '#', 'label']
            for j in range(int(nbondtypes)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l==0:
                        bondcoeffid = int(el[0])
                        bond[bondcoeffid] = {}
                        bond[bondcoeffid][str(k)] = {}
                    bond[bondcoeffid][str(k)] = (el[l+1])
                    #bond[bondcoeffid]['k'] = float(bond[bondcoeffid]['k'])
                    #bond[bondcoeffid]['r'] = float(bond[bondcoeffid]['r'])
            print "Bond Coeffs"
            print bond

        if "Angle Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns = ['k', 'theta', '#', 'label']
            for j in range(int(nangletypes)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        anglecoeffid = int(el[0])
                        angle[anglecoeffid] = {}
                        angle[anglecoeffid][str(k)] = {}
                    angle[anglecoeffid][str(k)] = (el[l+1])
            print "Angle Coeffs"
            print angle

        if "Dihedral Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns = ['a', 'b', 'c', 'd', '#', 'label']
            for j in range(int(ndhtypes)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        dhcoeffid = int(el[0])
                        dh[dhcoeffid] = {}
                        dh[dhcoeffid][str(k)] = {}
                    dh[dhcoeffid][str(k)] = (el[l+1])
            print "Dihedral Coeffs"
            print dh

        if "Improper Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns = ['a', 'b', 'c', 'd', '#', 'label']
            for j in range(int(nimptypes)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        impcoeffid = int(el[0])
                        imp[impcoeffid] = {}
                        imp[impcoeffid][str(k)] = {}
                    imp[impcoeffid][str(k)] = (el[l+1])
            print "Improper Coeffs"
            print imp


        if "Masses" in lines[i]:
            junk = lines[i+1]
            masscol = ['mass', 'type']
            for j in range(int(natomtypes)):
                el = lines[i+2+j].split()
                print el
                mass[el[3]] = float(el[1])
                #print mass[el[4]]

        if "Atoms" in lines[i]:
            junk = lines[i+1]
            columns = ['mol', 'type', 'charge', 'x', 'y', 'z', '#', 'atomname']
            print columns
            print columns[1]
            for j in range(int(natoms)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        atomid = int(el[0])
                        atomdata[atomid]= {}
                        atomdata[atomid][str(k)] = {}
                    atomdata[atomid][str(k)] = el[l+1]
                atomdata[atomid]['mol'] = int(atomdata[atomid]['mol'])
                atomdata[atomid]['type'] = int(atomdata[atomid]['type'])
                atomdata[atomid]['charge'] = float(atomdata[atomid]['charge'])
                atomdata[atomid]['x'] = float(atomdata[atomid]['x'])
                atomdata[atomid]['y'] = float(atomdata[atomid]['y'])
                atomdata[atomid]['z'] = float(atomdata[atomid]['z'])
                    #print atomnames
                """
                atomdata[atomid]['type'] = atomnames[atomdata[atomid]['type']]
                atomdata[atomid]['mol'] = atomnames[atomdata[atomid]['mol']]
                atomdata[atomid]['charge'] = atomnames[atomdata[atomid]['charge']]
                atomdata[atomid]['x'] = atomnames[atomdata[atomid]['x']]
                atomdata[atomid]['y'] = atomnames[atomdata[atomid]['y']]
                atomdata[atomid]['z'] = atomnames[atomdata[atomid]['z']]"""

        if "Bonds" in lines[i]:
            junk = lines[i+1]
            columns = ['bondtype', 'at1', 'at2']
            for j in range(int(nbonds)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        bondid = int(el[0])
                        bonddata[bondid]= {}
                        bonddata[bondid][str(k)] = {}
                    bonddata[bondid][str(k)] = int(el[l+1])
                    #print atomnames

        if "Angles" in lines[i]:
            junk = lines[i+1]
            columns = ['angletype', 'at1', 'at2', 'at3']
            for j in range(int(nangles)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        angleid = int(el[0])
                        angledata[angleid]= {}
                        angledata[angleid][str(k)] = {}
                    angledata[angleid][str(k)] = int(el[l+1])
                    #print atomnames

        if "Dihedrals" in lines[i]:
            junk = lines[i+1]
            columns = ['dhtype', 'at1', 'at2', 'at3', 'at4']
            for j in range(int(ndh)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        dhid = int(el[0])
                        dhdata[dhid]= {}
                        dhdata[dhid][str(k)] = {}
                    dhdata[dhid][str(k)] = int(el[l+1])
                    #print atomnames

        if "Impropers" in lines[i]:
            junk = lines[i+1]
            columns = ['imptype', 'at1', 'at2', 'at3', 'at4']
            for j in range(int(nimp)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        impid = int(el[0])
                        impdata[impid]= {}
                        impdata[impid][str(k)] = {}
                    impdata[impid][str(k)] = int(el[l+1])
                    #print atomnames

    return atomdata, bonddata, angledata, dhdata, impdata, mass, boxsize, numdata, pair, bond, angle, dh, imp

###########################################################
###########################################################
##  READ DATA FILE WITH NO FORCE FIELD PARAMETERS


def readpreFFData(fn):
    lines = readAll(fn)
    junk = ""
    atomdata = {}
    bonddata = {}
    angledata = {}
    dhdata = {}
    impdata = {}
    mass = {}
    boxsize = {}
    numdata = {}
    pair = {}
    bond = {}
    angle = {}
    dh = {}
    imp = {}

    numdata['natoms'] = {}
    numdata['nbonds'] = {}
    numdata['nangles'] = {}
    numdata['ndh'] = {}
    numdata['nimp'] = {}

    numdata['natomtypes'] = {}
    numdata['nbondtypes'] = {}
    numdata['nangletypes'] = {}
    numdata['ndhtypes'] = {}
    numdata['nimptypes'] = {}

    for i, x in enumerate(lines):
        if "atoms" in lines[i]:
            natoms = int(((lines[i].split("atoms")[0]).replace(" ","")).strip("\n"))
            numdata['natoms'] = natoms
        if "bonds" in lines[i]:
            nbonds = int(((lines[i].split("bonds")[0]).replace(" ","")).strip("\n"))
            numdata['nbonds'] = nbonds
        if "angles" in lines[i]:
            nangles = int(((lines[i].split("angles")[0]).replace(" ","")).strip("\n"))
            numdata['nangles'] = nangles
        if "dihedrals" in lines[i]:
            ndh = int(((lines[i].split("dihedrals")[0]).replace(" ","")).strip("\n"))
            numdata['ndh'] = ndh
        if "impropers" in lines[i]:
            nimp = int(((lines[i].split("impropers")[0]).replace(" ","")).strip("\n"))
            numdata['nimp'] = nimp
        if "atom types" in lines[i]:
            natomtypes = int(((lines[i].split("atom")[0]).replace(" ","")).strip("\n"))
            numdata['natomtypes'] = natomtypes
            print "natomtypes "+str(natomtypes)
        if "bond types" in lines[i]:
            nbondtypes = int(((lines[i].split("bond")[0]).replace(" ","")).strip("\n"))
            numdata['nbondtypes'] = nbondtypes
            print "nbondtypes "+str(nbondtypes)
        if "angle types" in lines[i]:
            nangletypes = int(((lines[i].split("angle")[0]).replace(" ","")).strip("\n"))
            numdata['nangletypes'] = nangletypes
            print "nangletypes "+str(nangletypes)
        if "dihedral types" in lines[i]:
            ndhtypes = int(((lines[i].split("dihedral")[0]).replace(" ","")).strip("\n"))
            numdata['ndhtypes'] = ndhtypes
            print "ndhtypes "+str(ndhtypes)
        if "improper types" in lines[i]:
            nimptypes = int(((lines[i].split("improper")[0]).replace(" ","")).strip("\n"))
            numdata['nimptypes'] = nimptypes
            print "nimptypes "+str(nimptypes)
        if "xlo xhi" in lines[i]:
            ln = lines[i].split()
            xlo = float(ln[0])
            xhi = float(ln[1])
            boxsize['xlo'] = {}
            boxsize['xhi'] = {}
            boxsize['xlo'] = xlo
            boxsize['xhi'] = xhi
        if "ylo yhi" in lines[i]:
            ln = lines[i].split()
            ylo = float(ln[0])
            yhi = float(ln[1])
            boxsize['ylo'] = {}
            boxsize['yhi'] = {}
            boxsize['ylo'] = ylo
            boxsize['yhi'] = yhi
        if "zlo zhi" in lines[i]:
            ln = lines[i].split()
            zlo = float(ln[0])
            zhi = float(ln[1])
            boxsize['zlo'] = {}
            boxsize['zhi'] = {}
            boxsize['zlo'] = zlo
            boxsize['zhi'] = zhi

        if "Pair Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ '#', 'num', 'label']
            for j in range(int(natomtypes)):
                el = lines[i+2+j].split()
                paircoeffid = int(el[1])
                pair[paircoeffid] = {}
                pair[paircoeffid]['label'] = {}
                pair[paircoeffid]['label'] = el[2]

                    #pair[paircoeffid]['eps'] = float(pair[paircoeffid]['eps'])
                    #pair[paircoeffid]['sig'] = float(pair[paircoeffid]['sig'])
            print "Pair Coeffs:"
            print pair

        if "Bond Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ '#', 'num', 'label']
            for j in range(int(nbondtypes)):
                el = lines[i+2+j].split()
                bondcoeffid = int(el[1])
                bond[bondcoeffid] = {}
                bond[bondcoeffid]['label'] = {}
                bond[bondcoeffid]['label'] = el[2]

                    #bond[bondcoeffid]['k'] = float(bond[bondcoeffid]['k'])
                    #bond[bondcoeffid]['r'] = float(bond[bondcoeffid]['r'])
            print "Bond Coeffs"
            print bond

        if "Angle Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ '#', 'num', 'label']
            for j in range(int(nangletypes)):
                el = lines[i+2+j].split()
                anglecoeffid = int(el[1])
                angle[anglecoeffid] = {}
                angle[anglecoeffid]['label'] = {}
                angle[anglecoeffid]['label'] = el[2]

            print "Angle Coeffs"
            print angle

        if "Dihedral Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ '#', 'num', 'label']
            for j in range(int(ndhtypes)):
                el = lines[i+2+j].split()
                dhcoeffid = int(el[1])
                dh[dhcoeffid] = {}
                dh[dhcoeffid]['label'] = {}
                dh[dhcoeffid]['label'] = el[2]

            print "Dihedral Coeffs"
            print dh

        if "Improper Coeffs" in lines[i]:
            junk = lines[i+1] # Reads blank line
            columns  = [ '#', 'num', 'label']
            for j in range(int(nimptypes)):
                el = lines[i+2+j].split()
                impcoeffid = int(el[1])
                imp[impcoeffid] = {}
                imp[impcoeffid]['label'] = {}
                imp[impcoeffid]['label'] = el[2]

            print "Improper Coeffs"
            print imp

        if "Atoms" in lines[i]:
            junk = lines[i+1]
            columns = ['mol', 'type', 'charge', 'x', 'y', 'z', '#', 'atomname']
            print columns
            print columns[1]
            for j in range(int(natoms)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        atomid = int(el[0])
                        atomdata[atomid]= {}
                        atomdata[atomid][str(k)] = {}
                    atomdata[atomid][str(k)] = el[l+1]
                atomdata[atomid]['mol'] = int(atomdata[atomid]['mol'])
                atomdata[atomid]['type'] = int(atomdata[atomid]['type'])
                atomdata[atomid]['charge'] = float(atomdata[atomid]['charge'])
                atomdata[atomid]['x'] = float(atomdata[atomid]['x'])
                atomdata[atomid]['y'] = float(atomdata[atomid]['y'])
                atomdata[atomid]['z'] = float(atomdata[atomid]['z'])
                    #print atomnames
                """
                atomdata[atomid]['type'] = atomnames[atomdata[atomid]['type']]
                atomdata[atomid]['mol'] = atomnames[atomdata[atomid]['mol']]
                atomdata[atomid]['charge'] = atomnames[atomdata[atomid]['charge']]
                atomdata[atomid]['x'] = atomnames[atomdata[atomid]['x']]
                atomdata[atomid]['y'] = atomnames[atomdata[atomid]['y']]
                atomdata[atomid]['z'] = atomnames[atomdata[atomid]['z']]"""

        if "Bonds" in lines[i]:
            junk = lines[i+1]
            columns = ['bondtype', 'at1', 'at2']
            for j in range(int(nbonds)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        bondid = int(el[0])
                        bonddata[bondid]= {}
                        bonddata[bondid][str(k)] = {}
                    bonddata[bondid][str(k)] = int(el[l+1])
                    #print atomnames

        if "Angles" in lines[i]:
            junk = lines[i+1]
            columns = ['angletype', 'at1', 'at2', 'at3']
            for j in range(int(nangles)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        angleid = int(el[0])
                        angledata[angleid]= {}
                        angledata[angleid][str(k)] = {}
                    angledata[angleid][str(k)] = int(el[l+1])
                    #print atomnames

        if "Dihedrals" in lines[i]:
            junk = lines[i+1]
            columns = ['dhtype', 'at1', 'at2', 'at3', 'at4']
            for j in range(int(ndh)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        dhid = int(el[0])
                        dhdata[dhid]= {}
                        dhdata[dhid][str(k)] = {}
                    dhdata[dhid][str(k)] = int(el[l+1])
                    #print atomnames

        if "Impropers" in lines[i]:
            junk = lines[i+1]
            columns = ['imptype', 'at1', 'at2', 'at3', 'at4']
            for j in range(int(nimp)):
                el = lines[i+2+j].split()
                for l,k in enumerate(columns):
                    if l == 0:
                        impid = int(el[0])
                        impdata[impid]= {}
                        impdata[impid][str(k)] = {}
                    impdata[impid][str(k)] = int(el[l+1])
                    #print atomnames

    return atomdata, bonddata, angledata, dhdata, impdata, mass, boxsize, numdata, pair, bond, angle, dh, imp



##########################################################
##########################################################



#########################################################
"""
getAt
last modified: 13/03/14
"""
def getAt(traj, atomdata, at1, at2):
    if at1 == 'first':
        at1  = sorted(atomdata.keys())[0]
    if at2 == 'last':
        at2 = sorted(atomdata.keys())[len(atomdata)-1]
    sortedatomnum = sorted(atomdata.keys())
    atrange = []
    for i,x in enumerate(sortedatomnum):
        if int(x) >= int(at1):
            if int(x) <= int(at2):
                atrange.append(x)
    return atrange
#######################################################

"""
getMol
last modified: 13/03/14
"""
def getMol(traj, atomdata, *mol):
    typerange = []
    sortedatomnum = sorted(atomdata.keys())
    print mol
    for i, x in enumerate(mol):
        for j, y in enumerate(sortedatomnum):
            pass
            if atomdata[y]["type"] == x:
                typerange.append(y)
    return typerange

"""
getTotMass
#get total mass of system
last modified: 13/03/14
"""
def getTotMass(traj, masses, natoms, tsrange):
    ts1 = sorted(traj.keys())[0] #first timestep
    totmass = 0
    atomnumrange = traj[ts1]['atom'].keys()
    for i in atomnumrange:
        totmass = totmass + masses[traj[ts1]['atom'][i]['type']]
    #for i in range(int(natoms)):
    #    totmass = totmass + masses[traj[ts1]['atom'][i+1]['type']]
    return totmass


################################################################

"""
getCOMTS
#calculate centre of mass
# Reads in timestep by timesteps
# returns the COM dictionary where the format is
#  COM[t,'xcom, ycom, zcom']
# !!! where t is the timestep from 1 to N, NOT the
# physical timestep number in simulation !!!
last modified: 13/03/14
"""
def getCOMts(lines, masses, natoms, header, atomnames):
    print 'Centre of Mass - start'
    tsrange = getTSrange(lines)
    COM = {}
    for i,t in enumerate(tsrange):
        datats = readTS(lines, header, natoms, t , atomnames)
        #print datats
        COM[t] = {}
        xf = 0.0
        yf = 0.0
        zf = 0.0
        if (i==0):
            atomnumrange = datats[tsrange[t]]['atom'].keys()
        #for j in range(int(natoms)):
        for j in atomnumrange:
            x = float(datats[tsrange[t]]['atom'][j]['x']) * float(masses[datats[tsrange[t]]['atom'][j]['type']])
            xf = xf + x
            y = float(datats[tsrange[t]]['atom'][j]['y']) * float(masses[datats[tsrange[t]]['atom'][j]['type']])
            yf = yf + y
            z = float(datats[tsrange[t]]['atom'][j]['z']) * float(masses[datats[tsrange[t]]['atom'][j]['type']])
            zf = zf + z
        if (i==0):
            totmass = getTotMass(datats, masses,natoms, tsrange)
        xcom = xf/totmass
        ycom = yf/totmass
        zcom = zf/totmass
        COM[t]['xcom'] = {}
        COM[t]['ycom'] = {}
        COM[t]['zcom'] = {}
        COM[t]['xcom'] = xcom
        COM[t]['ycom'] = ycom
        COM[t]['zcom'] = zcom
       # if
        if (t % 20 == 0):
            print 'COM step '+str(t)
    print 'Centre of Mass - complt'
    return COM
#############################################

"""
readcombmasses
#########################################
## READ COMBINED MASSES FOR DENSITY PROF
#########################################
last modified: 13/03/14
"""
def readcombmasses(fn, atommasses):
    lines=readAll(fn)
    combmass = {}
    for i,x in enumerate(lines):
        elems = lines[i].split()
        combmass[elems[0]] = {}
        combmass[elems[0]] = 0.0
        for j, y in enumerate(elems):
            combmass[elems[0]] += float(atommasses[y])
    return combmass


##########################################


####################################################################################
####################################################################################
##################### DEPRECATED VERSIONs ##########################################
####################################################################################

def readTraj(lines, header, natoms,tstot,atomnames):
    traj = {}
    j = 0
    k = 0
    tscount = 0
    for i, x in enumerate(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1].strip("\n"))
            traj[ts] = {}
            traj[ts]["atom"] = {}
            tscount +=1
            if (int(tscount) % 1 == 0):
                print "Read Timestep: " + str(tscount)+" / "+str(tstot)
        if "BOX BOUNDS pp pp pp" in lines[i]:
            boxx = map(float,lines[i+1].split())
            boxy = map(float, lines[i+2].split())
            boxz = map(float, lines[i+3].split())
            header = 9
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            traj[ts]["boxsize"] = (lx, ly, lz)
        if "xy xz yz pp pp pp" in lines[i]:
            boxx = map(float,lines[i+1].split())
            boxy = map(float, lines[i+2].split())
            boxz = map(float, lines[i+3].split())
            header = 9
            lx = boxx[1] - boxx[0]
            ly = boxy[1] - boxy[0]
            lz = boxz[1] - boxz[0]
            traj[ts]["boxsize"] = (lx, ly, lz)
        if lines[i].startswith("ITEM: ATOMS"):
            elems = lines[i].split()
            for j in range(natoms):
                el = lines[i+j+1].split()
                for l, k in enumerate(elems[2:]):
                    if l==0:
                        atomid=int(el[l])
                        traj[ts]["atom"][atomid]={}
                        continue
                    traj[ts]["atom"][atomid][k] = el[l]
                #Change from numerical atom type (from dump) to atom type name (from data file)
                traj[ts]["atom"][atomid]["type"] = atomnames[traj[ts]["atom"][atomid]["type"]]
    return traj

###########################################
## OBSOLETE - READS FROM WHOLE TRAJECTORY AT ONCE
def getCOM(traj, masses,natoms, tsrange):
    totmass = getTotMass(traj, masses,natoms)
    COM = {}
    tsrange = getTSrange()
    for i, t in enumerate(tsrange):
        COM[t] = {}
        xf = 0.0
        yf = 0.0
        zf = 0.0
        for j in range(int(natoms)):
            x = float(traj[t]['atom'][j+1]['x']) * float(masses[traj[t]['atom'][j+1]['type']])
            xf = xf + x
            y = float(traj[t]['atom'][j+1]['y']) * float(masses[traj[t]['atom'][j+1]['type']])
            yf = yf + y
            z = float(traj[t]['atom'][j+1]['z']) * float(masses[traj[t]['atom'][j+1]['type']])
            zf = zf + z
        xcom = xf/totmass
        ycom = yf/totmass
        zcom = zf/totmass
        COM[t]['xcom'] = {}
        COM[t]['ycom'] = {}
        COM[t]['zcom'] = {}
        COM[t]['xcom'] = xcom
        COM[t]['ycom'] = ycom
        COM[t]['zcom'] = zcom
    return COM
#########################################

def getVelProfTS(lines, masses, natoms, header, atomnames):
    print 'Calculate velocity profile'
    tsrange = getTSrange(lines)
    velprof = {}
    for i,t in enumerate(tsrange):
        datats = readTS(lines, header, natoms, t , atomnames)
        velprof[t] = {}
        vxf = 0.0
        vyf = 0.0
        vzf = 0.0
        if (i==0):
            atomnumrange = datats[tsrange[t]]['atom'].keys()
        #for j in range(int(natoms)):

           #NEEDS TO BE BINNED AT THIS STAGE - BASED ON z distance
        for j in atomnumrange:
            vx = float(datats[tsrange[t]]['atom'][j]['vx'])
            vxf += vx
            vy = float(datats[tsrange[t]]['atom'][j]['vy'])
            vyf += vy
            vz = float(datats[tsrange[t]]['atom'][j]['vz'])
            vzf += vz

        #if (i==0):
        #    totmass = getTotMass(datats, masses,natoms, tsrange)
        #xcom = xf/totmass
        #ycom = yf/totmass
        #zcom = zf/totmass
        velprof[t]['vx'] = {}
        velprof[t]['vy'] = {}
        velprof[t]['vz'] = {}
        velprof[t]['vz'] = vxf/float(len(atomrange))
        velprof[t]['vy'] = vyf/float(len(atomrange))
        velprof[t]['vz'] = vzf/float(len(atomrange))
        if (t % 20 == 0):
            print 'Vel. profile step '+str(t)
        print 'Vel. profile - complt'


