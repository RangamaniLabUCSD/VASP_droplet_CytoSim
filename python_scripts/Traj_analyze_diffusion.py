
nr = 4
######## Imports ########
#%matplotlib qt
import math
import numpy as np
import os
import pandas as pd
def generateemptylist(nrows,ncols):
    listret = []
    for r in range(nrows):
        lrow= []
        for c in range(ncols):
            lrow.append([])
        listret.append(lrow)
    return listret

import subprocess
class replicate:
    ridx = 0;
    Nsnaps = 0;
    snap = [];
    timevector = []
    def __init__(self, rid):
        self.ridx = rid
        self.Nsnaps = 0
        self.snap = [];
        self.timevector =[];
    def addsnapshot(self):
        self.snap.append(snapshot(1))


class snapshot:
    sidx = 0;
    filcoord=[]
    crossboundcoord=[];
    crossboundblobid = [];
    crossboundfilid = [];
    solidcoord = [];
    def __init__(self, sid):
        self.sidx = sid
        self.filcoord = [];
        self.crossboundcoord = [];
        self.crossboundblobid = [];
        self.crossboundfilid = [];
        self.solidcoord = [];

def readcytosimtraj (*args):
    if(len(args)):
        filepath = args[0]
        filename = 'objects.cmo'
    else:
        filepath = '/Users/aravind/Research/PostDoc/Research/cytosimfiles/Fgrow_3_2/kon_1e-1_koff_1e+0/'
        filename = 'objects.cmo'

    # Get number of solids
    test = subprocess.run(["./get_nsolids.sh","-f", filepath],
                        capture_output=True)
    nsolids = int(test.stdout)
    print("Nsolids in trajectory = "+str(nsolids), flush=True)
    r=[];
    Nruns =1 
    for  i in range(0,Nruns):
        r.append(replicate(i))
    printstatus = False;

    #filepath = '/Users/aravind/Research/PostDoc/Research/cytosimresults/sample_yossi/'
    #filename = 'objects.cmo'
    ridx = 0;
    sidx=-1;
    # Open file
    fptr = open(filepath+filename,'r')
    for line in fptr:
        if(printstatus):
            print(line)

        if('time' in line):
            t = float((line.split(' '))[1])
            r[ridx].timevector.append(t)
            r[ridx].addsnapshot();
            sidx = sidx + 1
            r[ridx].snap[sidx].solidcoord = [None]*nsolids

        elif('#section fiber' in line):
            fcoord=[];
            line = fptr.readline()
            if(printstatus):
                print(line)
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                if('f1' in line):
                    if(fcoord):
                        r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                        fcoord=[];
                elif(' ' in line[0]):
                    line = line.strip()
                    cstring = line.split(' ')
                    fcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                line = fptr.readline()
            #when it exists, if fcoord has not been recorded, record it.
            if(fcoord):
                r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                fcoord=[];
        # crosslinker hands
        if('#section solid' in line):
            line = fptr.readline()
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                #Get solid ID
                line = line.strip()
                cstring = line.split(' ')
                cstring = cstring[0].split(':')
                solidid = int(cstring[1])
                solidcoord = [];
                line = fptr.readline()
                while(line[0]!='d' and not('section' in line)):
                    #Get coord
                    line = line.strip()
                    cstring = line.split(' ')
                    solidcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                    line = fptr.readline()
                if(solidcoord):
                    r[ridx].snap[sidx].solidcoord[solidid-1] = np.array(solidcoord[0][:])
            
            # Bound crosslinker hands
            if('#section single A' in line):
                blobid = []
                line = fptr.readline()
                if(printstatus):
                    print(line)
                while(not('section' in line)):
                    if(printstatus):
                        print(line)
                    line = line.strip()
                    cstring = line.split(' ')
                    r[ridx].snap[sidx].crossboundfilid.append(int(cstring[1][1:len(cstring[3])]))
                    r[ridx].snap[sidx].crossboundblobid.append(int(cstring[3][1:len(cstring[3])]))
                    line = fptr.readline()
    fptr.close()
    print('Number of snapshots='+str(len(r[ridx].snap)))
    fptr.close()
    return r

def calculateMSD(solidcoordmat, MSD, MSD_snapcounter, MAXLAG):
    Nsnaps = len(solidcoordmat)
    for SREF in range(0, Nsnaps-1):
        coordref = solidcoordmat[SREF][:]
        maxstep = Nsnaps - SREF + 1
        for LREF in range(1,np.minimum(Nsnaps - SREF,MAXLAG)):
            coordlag = solidcoordmat[SREF+LREF][:]
            msd = np.mean(np.power(coordref - coordlag, 2))
            MSD[LREF] = MSD[LREF] + msd
            MSD_snapcounter[LREF] = MSD_snapcounter[LREF] + 1
    return [MSD, MSD_snapcounter]

def getsolidcoords(r, Nsnaps, Nsolids, minsnap):
    solidcoordmat = np.zeros((Nsnaps-minsnap,3*Nsolids))
    for SREF in range(minsnap,Nsnaps):
        solidcoordmat[SREF-minsnap,:]=np.reshape(r.snap[SREF].solidcoord,(1,3*Nsolids))
    return solidcoordmat

def execute():
    nreps = 3
    nr = 1
    nc = 6
    dirlist = ['N2k','N5k','N10k','N25k','N50k','N100k']
    minsnap = 1000
    MAXLAG = 950
    MSD_result = np.zeros((len(dirlist),MAXLAG))
    for cid in range(0,nc):
        MSD = np.zeros(MAXLAG)
        MSD_snapcounter = np.zeros(MAXLAG)
        for repid in range(0,nreps):
            dirname = dirlist[cid]+'/'+'r_'+str(repid)
            print(dirname,flush=True)
            r = readcytosimtraj('/scratch/achandrasekaran/cytosimfiles/test_diffusion/'+dirname+'/')
            NSNAPS = len(r[0].snap)
            NSOLIDS = len(r[0].snap[0].solidcoord)
            solidcoordmat = getsolidcoords(r[0], NSNAPS, NSOLIDS, minsnap)
            del r
            retval = calculateMSD(solidcoordmat, MSD, MSD_snapcounter, MAXLAG)
            MSD = retval[0]
            MSD_snapcounter = retval[1]
            MSD_snapcounter[0] = 1
        MSD_result[cid][:] = np.divide(MSD, MSD_snapcounter)
    pd.DataFrame(MSD_result).to_csv('sample.csv', index = dirlist)   

def getbinnedpwd(coord, bin_range, Nbinvec):
    Nbeads = len(coord)
    counter = 0
    pwd = np.zeros(int(Nbeads*(Nbeads -1)/2))
    for i in range(0,Nbeads):
        coordref = coord[i][:]
        for j in range(i+1,Nbeads):
            coordnext = coord[j][:]
            pwd[counter] = np.linalg.norm(coordnext-coordref)
            counter = counter + 1
    # Bin the distances
    bin_vector = []
    for Nbins in Nbinvec:
        X = 2*np.histogram(pwd, Nbins, bin_range)[0]
        bin_vector.append(X)
    return bin_vector
#  RDF
def executeRDF():
    Rdrop = 1
    Nedges = 101
    Nbins = Nedges - 1;
    nreps = 1
    nr = 1
    nc = 6
    dirlist = ['N2k','N5k','N10k','N25k','N50k','N100k']
    minsnap = 1000
    # Try 50, 100, 150, 200
    binNvec =np.array([50, 100, 150, 200])
    RDF_result_N = [];
    dR_N = []
    bin_edges_N = []
    bin_center_N = []
    bin_vector_S_N = []
    # set up bin edge vector
    # set up bin vector
    for bval in binNvec:
        RDF_result_N.append(np.zeros((len(dirlist),bval)))
        dR = 2*Rdrop/bval
        dR_N.append(dR)
        bedges = (np.linspace(0,2*Rdrop,bval+1));
        bin_edges_N.append(bedges)
        bin_center_N.append(bedges[0:bval:1]+dR/2)
        bin_vector_S_N.append(np.zeros(bval))
    for cid in range(0,nc):
        for repid in range(0,nreps):
            dirname = dirlist[cid]+'/'+'r_'+str(repid)
            print(dirname,flush=True)
            r = readcytosimtraj('/scratch/achandrasekaran/cytosimfiles/test_diffusion/'+dirname+'/')
            NSNAPS = len(r[0].snap)
            NSOLIDS = len(r[0].snap[0].solidcoord)
            snapcounter = 0
            # Get coordinates
            for SREF in range(minsnap,NSNAPS,50):
                coord = r[0].snap[SREF].solidcoord
                bin_range = np.array([0,2*Rdrop])
                bin_vector_S_list = getbinnedpwd(coord, bin_range, [50,100,150,200])
                for bREF in range(0,len(binNvec)):
                    bin_vector_S_N[bREF] = bin_vector_S_N[bREF] + np.cumsum(bin_vector_S_list[bREF])
                snapcounter = snapcounter + 1
        for bREF in range(0,len(binNvec)):
            Normalization_factor = np.power(bin_center_N[bREF],-2)*Rdrop**3/(3*snapcounter*NSOLIDS*NSOLIDS*dR_N[bREF])
            RDF_result_N[bREF][cid][:] = np.transpose(bin_vector_S_N[bREF]*Normalization_factor)
        for bREF in range(0,len(binNvec)):
            pd.DataFrame(RDF_result_N[bREF]).to_csv('RDF_'+dirlist[cid]+'_'+str(binNvec[bREF])+'.csv', index = dirlist)  
