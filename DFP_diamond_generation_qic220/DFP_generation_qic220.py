#!/usr/bin/python
########################################################################################################
######  Code is for generating DFP binary system with complemenary and same seqeunce length    #########
######  Modifiation is needed for self-complementary sequence and non-equal length sequence    #########
######  Interaction is excluded on same particle                                               #########
######  DNA strands is randomly grafted on core partilces                                      #########
######  Modification is needed to evenly graft AB strands                                      #########
######  Runfang Mao  5/25/2018    Lehigh University                                            #########  
########################################################################################################                               

#import MDAnalysis
import numpy as np
import math
import random
import copy
#import colloidsurface
from numpy import pi, cos, sin, arccos, arange

f = open("warm.data", "w")
setup = open("system_setup.dat", "r")

lattice_file = 'xyzfile/square.xyz'
numParticle_A = 0
numParticle_B = 0
DNAnum_A = 0
DNAnum_B = 0
percent_A = 0.0
percent_B = 0.0
totalnum = 0
box = 0.0
lz = 0.0
Rad_A = 0.0
Rad_B = 0.0
min_sep_A = 0.0
min_sep_B = 0.0
parsep = 50.0
dcolloid = 0.0
sequence_A = 'TT'
sequence_B = 'TT'
xyz = 'noprint'
dnaexist = 'yes'
print_lammps = 'no'

for line in setup.readlines():
    data = line.split()
    if type(vars()[data[0]]) == int:
        vars()[data[0]] = int(data[1])
    elif type(vars()[data[0]]) == float:
        vars()[data[0]] = float(data[1])
    elif type(vars()[data[0]]) == str:
        vars()[data[0]] = str(data[1])

numParticle = numParticle_A + numParticle_B

xtrialgroup = np.zeros(numParticle)
ytrialgroup = np.zeros(numParticle)
ztrialgroup = np.zeros(numParticle)

partcoordx = np.zeros([numParticle, totalnum])
partcoordy = np.zeros([numParticle, totalnum])
partcoordz = np.zeros([numParticle, totalnum])
partcoordx_B = np.zeros(totalnum)
partcoordy_B = np.zeros(totalnum)
partcoordz_B = np.zeros(totalnum)

Pos_A = []
Pos_B = []
Typeid_A = []
Typeid_B = []

Bondid = []
Bondid2 = []
Angleid = []


def cart2sph(x, y, z):
    XsqPlusYsq = x ** 2 + y ** 2
    r = math.sqrt(XsqPlusYsq + z ** 2)  # r
    elev = math.atan2(z, math.sqrt(XsqPlusYsq))  # theta
    az = math.atan2(y, x)  # phi
    return r, elev, az

def sph2cart(r, elev, az):
    x = r * math.cos(elev) * math.cos(az)
    y = r * math.cos(elev) * math.sin(az)
    z = r * math.sin(elev)
    return x, y, z

num_pts = DNAnum_A
indices = arange(0, num_pts, dtype=float) + 0.5

phi = arccos(1 - 2*indices/num_pts)
theta = pi * (1 + 5**0.5) * indices

cordx,cordy, cordz = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)

#for i in range(len(cordx)):
#    print 'check chekc chek',cordx[i],cordy[i],cordz[i], i

for i in range(num_pts):
        [r, elev, az] = cart2sph(cordx[i], cordy[i], cordz[i])
        r_A = r * Rad_A
        [x, y, z] = sph2cart(r_A, elev, az)
        partcoordx[0, i] = x
        partcoordy[0, i] = y
        partcoordz[0, i] = z
        r_B = r * Rad_B
        [x, y, z] = sph2cart(r_B, elev, az)
        partcoordx_B[i] = x
        partcoordy_B[i] = y
        partcoordz_B[i] = z

num_pts = totalnum - DNAnum_A
indices = arange(0, num_pts, dtype=float) + 0.5

phi = arccos(1 - 2*indices/num_pts)
theta = pi * (1 + 5**0.5) * indices

cordx,cordy, cordz = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)

for i in range(num_pts):
        [r, elev, az] = cart2sph(cordx[i], cordy[i], cordz[i])
        r_A = r * Rad_A
        [x, y, z] = sph2cart(r_A, elev, az)
        partcoordx[0, i+DNAnum_A] = x
        partcoordy[0, i+DNAnum_A] = y
        partcoordz[0, i+DNAnum_A] = z
        r_B = r * Rad_B
        [x, y, z] = sph2cart(r_B, elev, az)
        partcoordx_B[i+DNAnum_A] = x
        partcoordy_B[i+DNAnum_A] = y
        partcoordz_B[i+DNAnum_A] = z


C_count = 0
G_count = 0
seq_A1 = list(sequence_A)
lseq_A1 = len(sequence_A)
type_A1 = np.zeros(lseq_A1)
for i in range(0, lseq_A1):
    if seq_A1[i] == 'A':
        type_A1[i] = 3
    elif seq_A1[i] == 'T':
        type_A1[i] = 4
    elif seq_A1[i] == 'C':
        type_A1[i] = 5
        C_count += 1
    elif seq_A1[i] == 'G':
        type_A1[i] = 6
        G_count += 1

seq_B1 = list(sequence_B)
lseq_B1 = len(sequence_B)
type_B1 = np.zeros(lseq_B1)
for i in range(0, lseq_B1):
    if seq_B1[i] == 'A':
        type_B1[i] = 3
    elif seq_B1[i] == 'T':
        type_B1[i] = 4
    elif seq_B1[i] == 'C':
        type_B1[i] = 5
        C_count += 1
    elif seq_B1[i] == 'G':
        type_B1[i] = 6
        G_count += 1

if dnaexist == 'yes':
    ntypes = 8
    nbondtypes = 1
    nangletypes = 1
    ndihedraltypes = 0
else:
    ntypes = 1
    nbondtypes = 0
    nangletypes = 0
    ndihedraltypes = 0

### read DNA stacking type and hydrogen type ######
param_A1 = range(lseq_A1 - 1)
param_B1 = range(lseq_B1 - 1)

for i in range(1, lseq_A1):
    if (seq_A1[i - 1] == 'A' and seq_A1[i] == 'A'):
        param_A1[i - 1] = 'AA'
    elif (seq_A1[i - 1] == 'A' and seq_A1[i] == 'T') or (seq_A1[i - 1] == 'T' and seq_A1[i] == 'A'):
        param_A1[i - 1] = 'AT'
    elif (seq_A1[i - 1] == 'A' and seq_A1[i] == 'C') or (seq_A1[i - 1] == 'C' and seq_A1[i] == 'A'):
        param_A1[i - 1] = 'AC'
    elif (seq_A1[i - 1] == 'A' and seq_A1[i] == 'G') or (seq_A1[i - 1] == 'G' and seq_A1[i] == 'A'):
        param_A1[i - 1] = 'AG'
    elif (seq_A1[i - 1] == 'T' and seq_A1[i] == 'T'):
        param_A1[i - 1] = 'TT'
    elif (seq_A1[i - 1] == 'T' and seq_A1[i] == 'C') or (seq_A1[i - 1] == 'C' and seq_A1[i] == 'T'):
        param_A1[i - 1] = 'TC'
    elif (seq_A1[i - 1] == 'T' and seq_A1[i] == 'G') or (seq_A1[i - 1] == 'G' and seq_A1[i] == 'T'):
        param_A1[i - 1] = 'TG'
    elif (seq_A1[i - 1] == 'C' and seq_A1[i] == 'C'):
        param_A1[i - 1] = 'CC'
    elif (seq_A1[i - 1] == 'C' and seq_A1[i] == 'G') or (seq_A1[i - 1] == 'G' and seq_A1[i] == 'C'):
        param_A1[i - 1] = 'CG'
    elif (seq_A1[i - 1] == 'G' and seq_A1[i] == 'G'):
        param_A1[i - 1] = 'GG'

for i in range(1, lseq_B1):
    if (seq_B1[i - 1] == 'A' and seq_B1[i] == 'A'):
        param_B1[i - 1] = 'AA'
    elif (seq_B1[i - 1] == 'A' and seq_B1[i] == 'T') or (seq_B1[i - 1] == 'T' and seq_B1[i] == 'A'):
        param_B1[i - 1] = 'AT'
    elif (seq_B1[i - 1] == 'A' and seq_B1[i] == 'C') or (seq_B1[i - 1] == 'C' and seq_B1[i] == 'A'):
        param_B1[i - 1] = 'AC'
    elif (seq_B1[i - 1] == 'A' and seq_B1[i] == 'G') or (seq_B1[i - 1] == 'G' and seq_B1[i] == 'A'):
        param_B1[i - 1] = 'AG'
    elif (seq_B1[i - 1] == 'T' and seq_B1[i] == 'T'):
        param_B1[i - 1] = 'TT'
    elif (seq_B1[i - 1] == 'T' and seq_B1[i] == 'C') or (seq_B1[i - 1] == 'C' and seq_B1[i] == 'T'):
        param_B1[i - 1] = 'TC'
    elif (seq_B1[i - 1] == 'T' and seq_B1[i] == 'G') or (seq_B1[i - 1] == 'G' and seq_B1[i] == 'T'):
        param_B1[i - 1] = 'TG'
    elif (seq_B1[i - 1] == 'C' and seq_B1[i] == 'C'):
        param_B1[i - 1] = 'CC'
    elif (seq_B1[i - 1] == 'C' and seq_B1[i] == 'G') or (seq_B1[i - 1] == 'G' and seq_B1[i] == 'C'):
        param_B1[i - 1] = 'CG'
    elif (seq_B1[i - 1] == 'G' and seq_B1[i] == 'G'):
        param_B1[i - 1] = 'GG'

############################# LAMMPS data file ###################################################
pflag_A = int(percent_A * DNAnum_A)
pflag_B = int(percent_B * DNAnum_B)
ntotalatoms = numParticle * totalnum + numParticle_A * DNAnum_A * lseq_A1 * 2 + numParticle_B * DNAnum_B * lseq_B1 * 2
nbonds = numParticle_A * DNAnum_A * (lseq_A1 * 2 - 1) + numParticle_B * DNAnum_B * (
    lseq_B1 * 2 - 1) + numParticle_A * DNAnum_A + numParticle_B * DNAnum_B
nangles = numParticle_A * DNAnum_A * (lseq_A1 - 2) + numParticle_B * DNAnum_B * (lseq_B1 - 2)
# print 'percentage',percentage
# print 'pflag',pflag,'nNbase',nNbase,'nCbase',nCbase,'DNAnum',DNAnum

###############################  Atom coordinates  ##########################################################
keygroup_A = np.zeros(DNAnum_A)
keygroup_B = np.zeros(DNAnum_B)

#lattice_file = 'xyzfile/square.xyz'
core_x=[]
core_y=[]
core_z=[]
typelist = []
with open(lattice_file) as fh:
    next(fh)
    next(fh)
    for lines in fh:
        data = lines.split()
        typelist.append(data[0])
        core_x.append(float(data[1]))
        core_y.append(float(data[2]))
        core_z.append(float(data[3]))

#box = parsep*Lref
core_x = np.array(core_x)
core_y = np.array(core_y)
core_z = np.array(core_z)

xtrialgroup = core_x * parsep
ytrialgroup = core_y * parsep
ztrialgroup = core_z * parsep
xximage = np.around(xtrialgroup / box)
yyimage = np.around(ytrialgroup / box)
zzimage = np.around(ztrialgroup / box)
xtrialgroup -= box * xximage
ytrialgroup -= box * yyimage
ztrialgroup -= box * zzimage

file_surface = open('surface.xyz','w')
def atomprint(keygroup_A, keygroup_B):

    for j in range(totalnum):
        xxnew = partcoordx[0, j]  # surface particle
        yynew = partcoordy[0, j]
        zznew = partcoordz[0, j]
        xximage = np.around(xxnew / box)
        yyimage = np.around(yynew / box)
        zzimage = np.around(zznew / box)
        xxnew -= box * xximage
        yynew -= box * yyimage
        zznew -= box * zzimage
        Pos_A.append([xxnew, yynew, zznew])
        Pos_B.append([xxnew, yynew, zznew])

    ###################
    atomid = 3

    keygroup_A = np.arange(DNAnum_A)
    keygroup_A_s1 = np.array([0,31,10,14,12,21,8,25,29,4])
    keygroup_A_s2 = [x for x in keygroup_A if x not in keygroup_A_s1]

    for i in range(len(keygroup_A_s1)):
        keygroup_A[i] = keygroup_A_s1[i]
    for j in range(DNAnum_A - len(keygroup_A_s1)):
         keygroup_A[j + len(keygroup_A_s1)] = keygroup_A_s2[j]

    for n in range(DNAnum_A):
        gap =1.2
        diffx = partcoordx[0, int(keygroup_A[n])]
        diffy = partcoordy[0, int(keygroup_A[n])]
        diffz = partcoordz[0, int(keygroup_A[n])]
        diffrr = diffx * diffx + diffy * diffy + diffz * diffz
        diffr = math.sqrt(diffrr)
        norx = diffx / diffr  # claculate normal vector of picking particle
        nory = diffy / diffr
        norz = diffz / diffr
        for j in range(1, lseq_A1 * 2 + 1, 2):
            xxnew = partcoordx[0, int(keygroup_A[n])] + gap * norx
            yynew = partcoordy[0, int(keygroup_A[n])] + gap * nory
            zznew = partcoordz[0, int(keygroup_A[n])] + gap * norz
            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            Pos_A.append([xxnew, yynew, zznew])

            if diffx == 0:
                YY = 0.0
                XX = 1.0
            else:
                YY = 1 / (((diffy / diffx) ** 2 + 1) ** 0.5)
                XX = -YY * diffy / diffx

            xxnew = XX + partcoordx[0, int(keygroup_A[n])] + gap * norx
            yynew = YY + partcoordy[0, int(keygroup_A[n])] + gap * nory
            zznew = partcoordz[0, int(keygroup_A[n])] + gap * norz

            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            Pos_A.append([xxnew, yynew, zznew])
            atomid = atomid + 1

            gap = gap + 1.2

    keygroup_B = np.arange(DNAnum_B)

    count = 0
    for i in range(1, DNAnum_B):
        if i % int(DNAnum_B / pflag_B) == 0:
            count += 1

    keygroup_B_s1 = np.zeros(count)

    count = 0
    for i in range(1, DNAnum_B):
        if i % int(DNAnum_B / pflag_B) == 0:
            # print DNAnum_A/pflag_A, DNAnum_A,i,DNAnum_A%i
            keygroup_B_s1[count] = keygroup_B[i]
            count += 1
    keygroup_B_s2 = [x for x in keygroup_B if x not in keygroup_B_s1]

    for i in range(len(keygroup_B_s1)):
        keygroup_B[i] = keygroup_B_s1[i]
    for j in range(DNAnum_B - len(keygroup_B_s1)):
         keygroup_B[j + len(keygroup_B_s1)] = keygroup_B_s2[j]

    for n in range(DNAnum_B):
        gap =1.2
        diffx = partcoordx[0, int(keygroup_B[n])]
        diffy = partcoordy[0, int(keygroup_B[n])]
        diffz = partcoordz[0, int(keygroup_B[n])]
        diffrr = diffx * diffx + diffy * diffy + diffz * diffz
        diffr = math.sqrt(diffrr)
        norx = diffx / diffr  # claculate normal vector of picking particle
        nory = diffy / diffr
        norz = diffz / diffr
        for j in range(1, lseq_B1 * 2 + 1, 2):
            xxnew = partcoordx[0, int(keygroup_B[n])] + gap * norx
            yynew = partcoordy[0, int(keygroup_B[n])] + gap * nory
            zznew = partcoordz[0, int(keygroup_B[n])] + gap * norz
            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            Pos_B.append([xxnew, yynew, zznew])

            if diffx == 0:
                YY = 0.0
                XX = 1.0
            else:
                YY = 1 / (((diffy / diffx) ** 2 + 1) ** 0.5)
                XX = -YY * diffy / diffx

            xxnew = XX + partcoordx[0, int(keygroup_B[n])] + gap * norx
            yynew = YY + partcoordy[0, int(keygroup_B[n])] + gap * nory
            zznew = partcoordz[0, int(keygroup_B[n])] + gap * norz

            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            Pos_B.append([xxnew, yynew, zznew])
            atomid = atomid + 1

            gap = gap + 1.2

    for l in range(numParticle):
        print xtrialgroup[l], ytrialgroup[l], ztrialgroup[l]
        for i in range(totalnum):
            xxnew = Pos_A[i][0]+xtrialgroup[l]
            yynew = Pos_A[i][1]+ytrialgroup[l]
            zznew = Pos_A[i][2]+ztrialgroup[l]
            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            print xxnew,yynew,zznew

    for l in range(numParticle_A):
        for i in range(totalnum,len(Pos_A)):
            xxnew = Pos_A[i][0] + xtrialgroup[l]
            yynew = Pos_A[i][1] + ytrialgroup[l]
            zznew = Pos_A[i][2] + ztrialgroup[l]
            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            print xxnew, yynew, zznew

    for l in range(numParticle_B):
        for i in range(totalnum,len(Pos_B)):
            xxnew = Pos_B[i][0] + xtrialgroup[l+numParticle_A]
            yynew = Pos_B[i][1] + ytrialgroup[l+numParticle_A]
            zznew = Pos_B[i][2] + ztrialgroup[l+numParticle_A]
            xximage = np.around(xxnew / box)
            yyimage = np.around(yynew / box)
            zzimage = np.around(zznew / box)
            xxnew -= box * xximage
            yynew -= box * yyimage
            zznew -= box * zzimage
            print xxnew, yynew, zznew

    return keygroup_A, keygroup_B


###############################  type print  ##########################################################
def typeprint():
    for i in range(numParticle_A):
        print 'Acentral'
        #Typeid.append(11)
        for j in range(totalnum):
            print 1
            Typeid_A.append(1)
    for i in range(numParticle_B):
        print 'Bcentral'
        #Typeid.append(12)
        for j in range(totalnum):
            print 1
            Typeid_B.append(1)

    atomid = 3
    for l in range(numParticle_A):
        for j in range(pflag_A):
            for i in range(0, lseq_A1):
                print 2
                Typeid_A.append(2)
                print int(type_A1[i])
                Typeid_A.append(int(type_A1[i]))
                atomid += 1
        for j in range(DNAnum_A - pflag_A):
            for i in range(0, lseq_B1):
                print 2
                Typeid_A.append(2)
                print int(type_B1[i])
                Typeid_A.append(int(type_B1[i]))
                atomid += 1

    for l in range(numParticle_B):
        for j in range(pflag_B):
            for i in range(0, lseq_A1):
                print 2
                Typeid_B.append(2)
                print int(type_A1[i])
                Typeid_B.append(int(type_A1[i])+4)
                atomid += 1
        for j in range(DNAnum_B - pflag_B):
            for i in range(0, lseq_B1):
                print 2
                Typeid_B.append(2)
                print int(type_B1[i])
                Typeid_B.append(int(type_B1[i])+4)
                atomid += 1


###############################  body print  ##########################################################
def bodyprint():
    for l in range(numParticle):
        tagid = l * totalnum + l
        print tagid
        for j in range(totalnum):
            print tagid

    for l in range(numParticle_A):
        for j in range(pflag_A):
            for i in range(0, lseq_A1):
                print -1
                print -1
        for j in range(DNAnum_A - pflag_A):
            for i in range(0, lseq_B1):
                print -1
                print -1
    for l in range(numParticle_B):
        for j in range(pflag_B):
            for i in range(0, lseq_A1):
                print -1
                print -1
        for j in range(DNAnum_B - pflag_B):
            for i in range(0, lseq_B1):
                print -1
                print -1


###############################  mass print  ##########################################################
def massprint():
    for i in range(numParticle):
        print 1.0
        for j in range(totalnum):
            print 0.01
    for l in range(numParticle_A):
        for j in range(DNAnum_A):
            for i in range(1, lseq_A1 + 1):
                print 1.0
                print 1.0
    for l in range(numParticle_B):
        for j in range(DNAnum_B):
            for i in range(1, lseq_B1 + 1):
                print 1.0
                print 1.0


###############################  charge print  ##########################################################
def chargeprint():
    for i in range(numParticle):
        print 0.0
        for j in range(totalnum):
            print 0.0
    for l in range(numParticle_A):
        for j in range(DNAnum_A):
            for i in range(1, lseq_A1 + 1):
                print 0.0
                print 0.0
    for l in range(numParticle_B):
        for j in range(DNAnum_B):
            for i in range(1, lseq_B1 + 1):
                print 0.0
                print 0.0


################################ bonds ############################################################
def bondprint():
    t = 1
    for l in range(numParticle_A):
        for j in range(pflag_A):
            for i in range(0, lseq_A1 * 2, 2):
                bondini = numParticle * totalnum + l * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_A1 * 2
                if (i + 2 <= lseq_A1 * 2):
                    print 1, i + bondini + numParticle, i + 1 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 1 + bondini + numParticle])
                if (i + 3 <= lseq_A1 * 2):
                    print 1, i + bondini + numParticle, i + 2 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 2 + bondini + numParticle])
                if (i + 2 == lseq_A1 * 2):
                    t = t + 1
                else:
                    t = t + 2
        for j in range(DNAnum_A - pflag_A):
            for i in range(0, lseq_B1 * 2, 2):
                bondini = numParticle * totalnum + (l + 1) * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_B1 * 2
                if (i + 2 <= lseq_B1 * 2):
                    print 1, i + bondini + numParticle, i + 1 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 1 + bondini + numParticle])
                if (i + 3 <= lseq_B1 * 2):
                    print 1, i + bondini + numParticle, i + 2 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 2 + bondini + numParticle])
                if (i + 2 == lseq_B1 * 2):
                    t = t + 1
                else:
                    t = t + 2
    for l in range(numParticle_A):
        for i in range(DNAnum_A):
            bondini = numParticle * totalnum + l * pflag_A * lseq_A1 * 2 + l * (
            DNAnum_A - pflag_A) * lseq_B1 * 2 + i * lseq_A1 * 2 + 1
            print 1, '%d' % (keygroup_A[i] + l * totalnum + l + 1), bondini - 1 + numParticle
            Bondid2.append([keygroup_A[i] + l * totalnum + l + 1, bondini - 1 + numParticle])
            t = t + 1

    for l in range(numParticle_B):
        for j in range(pflag_B):
            for i in range(0, lseq_A1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + l * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_A1 * 2
                if (i + 2 <= lseq_A1 * 2):
                    print 1, i + bondini + numParticle, i + 1 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 1 + bondini + numParticle])
                if (i + 3 <= lseq_A1 * 2):
                    print 1, i + bondini + numParticle, i + 2 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 2 + bondini + numParticle])
                if (i + 2 == lseq_A1 * 2):
                    t = t + 1
                else:
                    t = t + 2
        for j in range(DNAnum_B - pflag_B):
            for i in range(0, lseq_B1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + (l + 1) * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_B1 * 2
                if (i + 2 <= lseq_B1 * 2):
                    print 1, i + bondini + numParticle, i + 1 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 1 + bondini + numParticle])
                if (i + 3 <= lseq_B1 * 2):
                    print 1, i + bondini + numParticle, i + 2 + bondini + numParticle
                    Bondid.append([i + bondini + numParticle, i + 2 + bondini + numParticle])
                if (i + 2 == lseq_B1 * 2):
                    t = t + 1
                else:
                    t = t + 2
    for l in range(numParticle_B):
        for i in range(DNAnum_B):
            bondini = numParticle * totalnum + numParticle_A * DNAnum_A * lseq_A1 * 2 + l * DNAnum_B * lseq_B1 * 2 + i * lseq_B1 * 2 + 1
            print 1, '%d' % (
                keygroup_B[i] + (l + numParticle_A) * totalnum + l + numParticle_A + 1), bondini - 1 + numParticle
            Bondid2.append(
                [keygroup_B[i] + (l + numParticle_A) * totalnum + l + numParticle_A + 1, bondini - 1 + numParticle])
            t = t + 1


############################# special pair ###########################
def specialstackprint():
    for l in range(numParticle_A):
        for j in range(pflag_A):
            count = 0
            for i in range(0, lseq_A1 * 2, 2):
                bondini = numParticle * totalnum + l * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_A1 * 2
                if (i + 3 <= lseq_A1 * 2):
                    print '%s' % param_A1[count], i + 1 + bondini + numParticle, i + 3 + bondini + numParticle
                    count += 1
        for j in range(DNAnum_A - pflag_A):
            count = 0
            for i in range(0, lseq_B1 * 2, 2):
                bondini = numParticle * totalnum + (l + 1) * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_B1 * 2
                if (i + 3 <= lseq_B1 * 2):
                    print '%s' % param_B1[count], i + 1 + bondini + numParticle, i + 3 + bondini + numParticle
                    count += 1

    for l in range(numParticle_B):
        for j in range(pflag_B):
            count = 0
            for i in range(0, lseq_A1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + l * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_A1 * 2
                if (i + 3 <= lseq_A1 * 2):
                    print '%s' % param_A1[count], i + 1 + bondini + numParticle, i + 3 + bondini + numParticle
                    count += 1
        for j in range(DNAnum_B - pflag_B):
            count = 0
            for i in range(0, lseq_B1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + (l + 1) * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_B1 * 2
                if (i + 3 <= lseq_B1 * 2):
                    print '%s' % param_B1[count], i + 1 + bondini + numParticle, i + 3 + bondini + numParticle
                    count += 1


######################
def specialwcaprint():
    for l in range(numParticle_A):
        for j in range(pflag_A):
            count = 0
            for i in range(0, lseq_A1 * 2, 2):
                bondini = numParticle * totalnum + l * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_A1 * 2
                if (i + 5 <= lseq_A1 * 2):
                    if (seq_A1[count] == 'A' and seq_A1[count + 2] == 'T') or (
                                    seq_A1[count] == 'T' and seq_A1[count + 2] == 'A'):
                        print 'hAT', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                    if (seq_A1[count] == 'C' and seq_A1[count + 2] == 'G') or (
                                    seq_A1[count] == 'G' and seq_A1[count + 2] == 'C'):
                        print 'hCG', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                count += 1
        for j in range(DNAnum_A - pflag_A):
            count = 0
            for i in range(0, lseq_B1 * 2, 2):
                bondini = numParticle * totalnum + (l + 1) * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_B1 * 2
                if (i + 5 <= lseq_B1 * 2):
                    if (seq_B1[count] == 'A' and seq_B1[count + 2] == 'T') or (
                                    seq_B1[count] == 'T' and seq_B1[count + 2] == 'A'):
                        print 'hAT', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                    if (seq_B1[count] == 'C' and seq_B1[count + 2] == 'G') or (
                                    seq_B1[count] == 'G' and seq_B1[count + 2] == 'C'):
                        print 'hCG', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                count += 1

    for l in range(numParticle_B):
        for j in range(pflag_B):
            count = 0
            for i in range(0, lseq_A1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + l * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_A1 * 2
                if (i + 5 <= lseq_A1 * 2):
                    if (seq_A1[count] == 'A' and seq_A1[count + 2] == 'T') or (
                                    seq_A1[count] == 'T' and seq_A1[count + 2] == 'A'):
                        print 'hAT', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                    if (seq_A1[count] == 'C' and seq_A1[count + 2] == 'G') or (
                                    seq_A1[count] == 'G' and seq_A1[count + 2] == 'C'):
                        print 'hCG', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                count += 1
        for j in range(DNAnum_B - pflag_B):
            count = 0
            for i in range(0, lseq_B1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + (l + 1) * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_B1 * 2
                if (seq_B1[count] == 'A' and seq_B1[count + 2] == 'T') or (
                                seq_B1[count] == 'T' and seq_B1[count + 2] == 'A'):
                    print 'hAT', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                if (seq_B1[count] == 'C' and seq_B1[count + 2] == 'G') or (
                                seq_B1[count] == 'G' and seq_B1[count + 2] == 'C'):
                    print 'hCG', i + 1 + bondini + numParticle, i + 5 + bondini + numParticle
                count += 1


######################
def specialwca_excluding_print():
    C_id = []
    G_id = []
    for l in range(numParticle_A):
        for j in range(pflag_A):
            count = 0
            for i in range(0, lseq_A1 * 2, 2):
                bondini = numParticle * totalnum + l * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_A1 * 2
                if i <= lseq_A1 * 2:
                    if (seq_A1[count] == 'C'):
                        # print 'hC', i + 1 + bondini + numParticle
                        C_id.append(i + 1 + bondini + numParticle)
                    if (seq_A1[count] == 'G'):
                        # print 'hG', i + 1 + bondini + numParticle
                        G_id.append(i + 1 + bondini + numParticle)
                count += 1
        for j in range(DNAnum_A - pflag_A):
            count = 0
            for i in range(0, lseq_B1 * 2, 2):
                bondini = numParticle * totalnum + (l + 1) * pflag_A * lseq_A1 * 2 + l * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2 + j * lseq_B1 * 2
                if i <= lseq_B1 * 2:
                    if (seq_B1[count] == 'C'):
                        # print 'hC', i + 1 + bondini + numParticle
                        C_id.append(i + 1 + bondini + numParticle)
                    if (seq_B1[count] == 'G'):
                        # print 'hG', i + 1 + bondini + numParticle
                        G_id.append(i + 1 + bondini + numParticle)
                count += 1

    C_position = np.reshape(C_id, (numParticle_A, pflag_A, C_count))
    G_position = np.reshape(G_id, (numParticle_A, DNAnum_A - pflag_A, G_count))
    # print C_position,G_position
    if C_count != 0 and G_count != 0:
        for i in range(numParticle_A):
            for j in range(pflag_A):
                for k in range(DNAnum_A - pflag_A):
                    for kk in range(C_count):
                        for kkk in range(G_count):
                            print 'hCG', C_position[i][j][kk], G_position[i][k][kkk]  # ,kk,kkk

    C_id = []
    G_id = []
    for l in range(numParticle_B):
        for j in range(pflag_B):
            count = 0
            for i in range(0, lseq_A1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + l * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_A1 * 2
                if i <= lseq_A1 * 2:
                    if (seq_A1[count] == 'C'):
                        # print 'hC', i + 1 + bondini + numParticle
                        C_id.append(i + 1 + bondini + numParticle)
                    if (seq_A1[count] == 'G'):
                        # print 'hG', i + 1 + bondini + numParticle
                        G_id.append(i + 1 + bondini + numParticle)
                count += 1
        for j in range(DNAnum_B - pflag_B):
            count = 0
            for i in range(0, lseq_B1 * 2, 2):
                particle_ini = numParticle * totalnum + numParticle_A * pflag_A * lseq_A1 * 2 + numParticle_A * (
                    DNAnum_A - pflag_A) * lseq_B1 * 2
                bondini = particle_ini + (l + 1) * pflag_B * lseq_A1 * 2 + l * (
                    DNAnum_B - pflag_B) * lseq_B1 * 2 + j * lseq_B1 * 2
                if i <= lseq_B1 * 2:
                    if (seq_B1[count] == 'C'):
                        # print 'hC', i + 1 + bondini + numParticle
                        C_id.append(i + 1 + bondini + numParticle)
                    if (seq_B1[count] == 'G'):
                        # print 'hG', i + 1 + bondini + numParticle
                        G_id.append(i + 1 + bondini + numParticle)
                count += 1
    C_position = np.reshape(C_id, (numParticle_B, pflag_B, C_count))
    G_position = np.reshape(G_id, (numParticle_B, DNAnum_B - pflag_B, G_count))
    # print C_position,G_position
    if C_count != 0 and G_count != 0:
        for i in range(numParticle_B):
            for j in range(pflag_B):
                for k in range(DNAnum_B - pflag_B):
                    for kk in range(C_count):
                        for kkk in range(G_count):
                            print 'hCG', C_position[i][j][kk], G_position[i][k][kkk]  # ,kk,kkk


############################### angle ####################################################################
def angleprint():
    t = 1
    for l in range(numParticle_A):
        #      print l, 'Particle'
        for j in range(DNAnum_A):
            #        print j, 'DNA'
            for i in range(0, lseq_A1 * 2, 2):
                bondini = j * lseq_A1 * 2 + l * DNAnum_A * lseq_A1 * 2 + numParticle * totalnum
                if (i + 5 <= lseq_A1 * 2):
                    print 1, i + bondini + numParticle, i + 2 + bondini + numParticle, i + 4 + bondini + numParticle
                    Angleid.append(
                        [i + bondini + numParticle, i + 2 + bondini + numParticle, i + 4 + bondini + numParticle])
                if (i + 5 == lseq_A1 * 2 - 1):
                    t = t - 1
                else:
                    t = t + 1

    for l in range(numParticle_B):
        #      print l, 'Particle'
        for j in range(DNAnum_B):
            #        print j, 'DNA'
            for i in range(0, lseq_B1 * 2, 2):
                bondini = j * lseq_B1 * 2 + numParticle_A * DNAnum_A * lseq_A1 * 2 + l * DNAnum_B * lseq_B1 * 2 + numParticle * totalnum
                if (i + 5 <= lseq_B1 * 2):
                    print 1, i + bondini + numParticle, i + 2 + bondini + numParticle, i + 4 + bondini + numParticle
                    Angleid.append(
                        [i + bondini + numParticle, i + 2 + bondini + numParticle, i + 4 + bondini + numParticle])
                if (i + 5 == lseq_B1 * 2 - 1):
                    t = t - 1
                else:
                    t = t + 1


############################ HOOMD input file ###################################################################


print '<?xml version ="1.0" encoding ="UTF-8" ?>'
print '<hoomd_xml version="1.4">'
print '<configuration time_step="0" natoms="%s">' % (ntotalatoms + numParticle)

print '<type num="%s">' % (ntotalatoms + numParticle)
typeprint()
print '</type>'

print '<body>'
bodyprint()
print '</body>'

print '<mass num="%s">' % (ntotalatoms + numParticle)
massprint()
print '</mass>'

# print '<charge num="%s">' %ntotalatoms
# chargeprint()
# print '</charge>'

print '<position num="%s">' % (ntotalatoms + numParticle)
keygroup_A, keygroup_B = atomprint(keygroup_A, keygroup_B)
print '</position>'

print '<bond num="%s">' % nbonds
bondprint()
print '</bond>'

print '<angle num="%s">' % nangles
angleprint()
print '</angle>'

print '<pair>'
specialstackprint()
specialwcaprint()
#specialwca_excluding_print()
print '</pair>'

print '<box lx="%s" ly="%s" lz="%s" />' % (box, box, lz)
print '</configuration>'
print '</hoomd_xml>'


#################### Lammps input file ##############################
def print_to_lammps():
    print >> f, 'LAMMPS DNA data file'
    print >> f, ''

    print >> f, ntotalatoms + numParticle, ' atoms'
    print >> f, nbonds, ' bonds'
    print >> f, nangles, ' angles'
    print >> f, 0, ' dihedrals'
    print >> f, 0, ' impropers'
    print >> f, ''

    print >> f, ntypes+4, ' atom types'
    print >> f, nbondtypes, ' bond types'
    print >> f, nangletypes, ' angle types'
    print >> f, ndihedraltypes, ' dihedral types'
    print >> f, 0, ' improper types'
    print >> f, ''

    print >> f, 0, box, ' xlo xhi'
    print >> f, 0, box, ' ylo yhi'
    print >> f, 0, box, ' zlo zhi'

    print >> f, ''
    print >> f, 'Masses'
    print >> f, ''

    print >> f, 1, 0.01

    for i in range(2, 10 + 1):
        print >> f, i, 1
    for b in range(0, numParticle):
        print >>f, 10+b+1, 0.01

    print >> f, ''
    print >> f, 'Atoms'
    print >> f, ''

    for i in range(totalnum):
        print >> f,i + 1, 1, Typeid_A[i], 0.0, Pos_A[i][0], Pos_A[i][1], Pos_A[i][2]
    for i in range(totalnum):
        print >> f, i+totalnum + 1, 2, Typeid_B[i], 0.0, Pos_B[i][0]+parsep, Pos_B[i][1]+parsep, Pos_B[i][2]+parsep
    for i in range(totalnum,len(Pos_A)):
        print >> f,i + totalnum + 1, 3, Typeid_A[i], 0.0, Pos_A[i][0], Pos_A[i][1], Pos_A[i][2]
    for i in range(totalnum , len(Pos_B)):
        print >> f, i + len(Pos_A)+1, 3, Typeid_B[i], 0.0, Pos_B[i][0]+parsep, Pos_B[i][1]+parsep, Pos_B[i][2]+parsep
    print >>f, ntotalatoms+numParticle-1,1,11,0.0,0.0,0.0,0.0
    print >>f, ntotalatoms+numParticle,2,12,0.0,0.0+parsep,0.0+parsep,0.0+parsep

    #print >> f, ''
    #print >> f, 'Bonds'
    #print >> f, ''
    #for i in range(len(Bondid)):
    #    print >> f, i + 1, 1, int(Bondid[i][0]) + 1, int(Bondid[i][1]) + 1
    #for i in range(len(Bondid2)):
    #    print >> f, i + 1 + len(Bondid), 1, int(Bondid2[i][0]) + 1, int(Bondid2[i][1]) + 1

    #print >> f, ''
    #print >> f, 'Angles'
    #print >> f, ''
    #for i in range(len(Angleid)):
    #    print >> f, i + 1, 1, Angleid[i][0] + 1, Angleid[i][1] + 1, Angleid[i][2] + 1


if 'yes' == print_lammps:
    print_to_lammps()
