######## Imports ########
#%matplotlib qt
####### This python scipt helps analyze the trajectories from CytoSim
import math
import numpy as np
import os
import pandas as pd
import icosphere as ic
from datetime import datetime
import os
def utilityN1e3(kbind,kunbind):
    kbindlist_set1 = [1e-3, 1e-2, 1e-1, 1e+0]
    kunbindlist_set1 = [1e-3, 1e-2, 1e-1, 1e+0]
    kbindlist_set2 = [1e-3, 1e-2, 1e-1, 1e+0]
    kunbindlist_set2 = [1e-5, 1e-4, 1e+1]
    kbindlist_set3 = [1e-4]
    kunbindlist_set3 = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0, 1e+1]
    kbindlist_set4 = [0.1]
    kunbindlist_set4 = [0.02, 0.04, 0.06, 0.08, 0.2, 0.4, 0.6, 0.8]
    kbindlist_set5 = [1.0]
    kunbindlist_set5 = [0.02, 0.04, 0.06, 0.08, 0.2, 0.4, 0.6, 0.8]
    kbindlist_set6 = [1e+1, 1e+2, 1e+3]
    kunbindlist_set6 = [1e-2, 1e-1, 1e+0, 1e+1]
    kbindlist_set7 = [1e+0]
    kunbindlist_set7 = [2.0,4.0,6.0,8.0]
    nr = -1
    if kbind in kbindlist_set1 and kunbind in kunbindlist_set1:
        nr = len(kbindlist_set1); nc = len(kunbindlist_set1)
        foldername = 'Fgrow_10_3_N_1e3_5reps'
        ridx = kbindlist_set1.index(kbind)
        cidx = kunbindlist_set1.index(kunbind)  
    elif kbind in kbindlist_set2 and kunbind in kunbindlist_set2:
        nr = len(kbindlist_set2); nc = len(kunbindlist_set2)
        foldername = 'Fgrow_10_3_N_1e3_5reps_set2'
        ridx = kbindlist_set2.index(kbind)
        cidx = kunbindlist_set2.index(kunbind)  
    elif kbind in kbindlist_set3 and kunbind in kunbindlist_set3:
        nr = len(kbindlist_set3); nc = len(kunbindlist_set3)
        foldername = 'Fgrow_10_3_N_1e3_5reps_set3'
        ridx = kbindlist_set3.index(kbind)
        cidx = kunbindlist_set3.index(kunbind)  
    elif kbind in kbindlist_set4 and kunbind in kunbindlist_set4:
        nr = len(kbindlist_set4); nc = len(kunbindlist_set4)
        foldername = 'Fgrow_10_3_N_1e3_5reps_set4'
        ridx = kbindlist_set4.index(kbind)
        cidx = kunbindlist_set4.index(kunbind)  
    elif kbind in kbindlist_set5 and kunbind in kunbindlist_set5:
        nr = len(kbindlist_set5); nc = len(kunbindlist_set5)
        foldername = 'Fgrow_10_3_N_1e3_5reps_set5'
        ridx = kbindlist_set5.index(kbind)
        cidx = kunbindlist_set5.index(kunbind)  
    elif kbind in kbindlist_set6 and kunbind in kunbindlist_set6:
        nr = len(kbindlist_set6); nc = len(kunbindlist_set6)
        foldername = 'Fgrow_10_3_N_1e3_5reps_set6'
        ridx = kbindlist_set6.index(kbind)
        cidx = kunbindlist_set6.index(kunbind)  
    elif kbind in kbindlist_set7 and kunbind in kunbindlist_set7:
        nr = len(kbindlist_set7); nc = len(kunbindlist_set7)
        foldername = 'Fgrow_10_3_N_1e3_5reps_set7'
        ridx = kbindlist_set7.index(kbind)
        cidx = kunbindlist_set7.index(kunbind)  
    if nr==-1:
        return []
    Rval = nc*ridx+cidx
    return [foldername,Rval]
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

def readcytosimreport (*args):
    tag = ''
    if(len(args)):
        filepath = args[0]
        filename = 'fiber_points.cmo'
        filename2 = 'single.cmo'
        if(len(args)==2):
            tag= args[1]
    else:
        print('ERROR!!!!')
        return 0;

    r=[];
    Nruns =1 
    for  i in range(0,Nruns):
        r.append(replicate(i))
    printstatus = False;
    recordStatus = False;
    ridx = 0;
    sidx=-1;
    fcoord=[];
    # Open file
    fptr = open(filepath+filename,'r')
    for line in fptr:
        if(printstatus):
            print(line)

        if('time' in line):
            t = float((line.split(' '))[2])
            r[ridx].timevector.append(t)
            r[ridx].addsnapshot();
            sidx = sidx + 1

        elif('fiber' in line):
            fcoord=[];
            line = fptr.readline()
            if(printstatus):
                print(line)
            while(not('end' in line)):
                if(printstatus):
                    print(line)
                if(fcoord and 'fiber' in line):
                    r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                    fcoord=[];
                if(' ' in line[0]):
                    line = line.strip()
                    cstring = line.split()
                    fcoord.append([float(cstring[1]), float(cstring[2]), float(cstring[3])])
                line = fptr.readline()
            #when it exists, if fcoord has not been recorded, record it.
            if(fcoord):
                r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                fcoord=[];
    fptr.close()
    print('Number of snapshots='+str(len(r[ridx].snap)))
    #### READ CROSSLINKERS
    sidx=-1;
    # Open file
    fptr = open(filepath+filename2,'r', encoding="utf8", errors="ignore")
    for line in fptr:
        if(printstatus):
            print(line)
        if('time' in line):
            sidx = sidx + 1
        if('class' in line):
            while(not('end' in line)):
                line= fptr.readline()
                if not('end' in line):
                    line = line.strip()
                    cstring = line.split()
                    filID = int(cstring[5])
                    if filID>0:
                        r[ridx].snap[sidx].crossboundfilid.append(filID)
                        r[ridx].snap[sidx].crossboundblobid.append(int(int(cstring[1])/4))
            
    fptr.close()
    return r

def readcytosimtraj (*args):
    tag = ''
    if(len(args)):
        filepath = args[0]
        filename = 'objects.cmo'
        if(len(args)==2):
            tag= args[1]
    else:
        filepath = '/Users/aravind/Research/PostDoc/Research/cytosimfiles/Fgrow_3_2/kon_1e-1_koff_1e+0/'
        filename = 'objects.cmo'
    
    if os.path.exists(filepath+'fiber_points.cmo'):
        print('Using cytosim report')
        r=readcytosimreport (*args);
        return r;
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
    recordStatus = False;
    #filepath = '/Users/aravind/Research/PostDoc/Research/cytosimresults/sample_yossi/'
    #filename = 'objects.cmo'
    ridx = 0;
    sidx=-1;
    # Open file
    fptr = open(filepath+filename,'r', encoding="utf8", errors="ignore")
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
                if('f2' in line):
                    recordStatus = False;
                if('f1' in line):
                    if(fcoord):
                        r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                        fcoord=[];
                    recordStatus = True;
                elif(' ' in line[0] and recordStatus):
                    line = line.strip()
                    cstring = line.split(' ')
                    fcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                line = fptr.readline()
            #when it exists, if fcoord has not been recorded, record it.
            if(fcoord):
                r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                fcoord=[];
        # crosslinker hands
        if('#section solid' in line and 'filonly' != tag):
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
            if('#section single A' in line and 'filonly' != tag):
                blobid = []
                line = fptr.readline()
                if(printstatus):
                    print(line)
                while(not('section' in line)):
                    if(printstatus):
                        print(line)
                    if line[0]=='w':
                        line = line.strip()
                        cstring = line.split(' ')
                        r[ridx].snap[sidx].crossboundfilid.append(int(cstring[1][1:len(cstring[1])]))
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

def getfillength(fc):
    Nbeads = (np.shape(fc))[0]
    L = 0
    for i in range(0,Nbeads-1):
        L = L + np.linalg.norm(fc[i]-fc[i+1])
    return L

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

def getfracoccupiedtimeseries(N, dirname, meshvar, Nlists, unitnormal, Ntriangles, r, deltasnap, outfilename):
    Nsnaps = len(r[0].snap)
    foccupied=[];
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    fvec = np.zeros((2,Ndatapoints))
    fvec[0][:]=np.arange(0,Nsnaps,deltasnap)
    for SREF in range(0,Nsnaps,deltasnap):
        print(SREF,flush=True)
        filcoord = r[0].snap[SREF].filcoord 
        actincounter = ic.generatedensityfield(meshvar, unitnormal, filcoord, False, ic.SearchAlgoType.LISTEDSEARCH, Nlists)
        #Order Parameter #3
        is_all_zero = np.argwhere((actincounter == 0.0))
        #foccupied.append(len(is_all_zero))
        fvec[1][SREF] = 1-len(is_all_zero)/Ntriangles
    pd.DataFrame(fvec, index=['Time','Frac_occupied']).to_csv(outfilename)

def getvaspprops(rset, deltasnap, N, Rval, outputfile):
    MEANID = 0
    STDID = 1
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    bar_min_snap = int(0.95*Nsnaps)
    bar_data_raw = np.array([])
    collectstatus = False
    datawritestatus = False
    datamatrix = np.zeros((19,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    crosslinker_n_fil_bar = [[],[],[],[]]
    # Nvasp = np.zeros((2,Ndatapoints))
    # Nvalency = np.zeros((2,Ndatapoints))
    # Nfil_per_solid = np.zeros((2,Ndatapoints))
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        nvasp_tmp = np.array([])
        nvalency_tmp = np.array([])
        Nfil_per_solid_tmp = np.array([])
        if SREF>=bar_min_snap:
            collectstatus = True
        countermat = np.zeros((4,len(rset)))
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                crossboundfilid = np.array(rset[runidx].snap[SREF].crossboundfilid)
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid
                if(len(crossboundblobid)):
                    datawritestatus = True
                unique, counts = np.unique(crossboundblobid, return_counts=True)
                nvasp_tmp = np.append(nvasp_tmp,len(unique))
                nvalency_tmp = np.append(nvalency_tmp,np.array(counts))
                # Get the filament IDs for each solid and calculate Nhands/Nfil
                # This helps you understand the number of hands in each solid that are bound to the same filament.
                # The distribution of Nhands/Nfil - closer to 1 suggests that each bound hand is bound to a different filament
                # >1 suggests that multiple hands of the solid are bound to the same fil ID
                # <1 is not possible.
                for uiter, ublobid in enumerate(unique):
                    #Find the locs using argwhere
                    locs = np.argwhere(crossboundblobid==ublobid)
                    #Corresponding filament IDs that share this 
                    #crosslinker in the bond
                    filid_solid = crossboundfilid[locs]
                    #Get filament IDs and the counts
                    unique_fid, counts_fil = np.unique(filid_solid, return_counts=True)
                    #Get the total number of unique filaments
                    lval = len(unique_fid)
                    countermat[lval-1][runidx] = countermat[lval-1][runidx] + 1
                    Nfil_per_solid_tmp = np.append(Nfil_per_solid_tmp,lval)
                    # Go through unique crosslinkers
                    # Find 
                if collectstatus:
                    bar_data_raw = np.append(bar_data_raw, Nfil_per_solid_tmp)
                    for pos in range(0,4):
                        crosslinker_n_fil_bar[pos].append(countermat[pos][runidx])

        meanvec = np.mean(countermat,1)
        sstdvec = np.std(countermat,1)


        if(datawritestatus):
            datamatrix[1][scounter]= np.mean(nvasp_tmp)
            datamatrix[2][scounter]= np.std(nvasp_tmp)
            datamatrix[3][scounter]= np.mean(nvalency_tmp)
            datamatrix[4][scounter]= np.std(nvalency_tmp)
            datamatrix[5][scounter]= np.mean(Nfil_per_solid_tmp)
            datamatrix[6][scounter]= np.std(Nfil_per_solid_tmp)
            datamatrix[9][scounter]= meanvec[0]
            datamatrix[10][scounter]= sstdvec[0]
            datamatrix[11][scounter]= meanvec[1]
            datamatrix[12][scounter]= sstdvec[1]
            datamatrix[13][scounter]= meanvec[2]
            datamatrix[14][scounter]= sstdvec[2]
            datamatrix[15][scounter]= meanvec[3]
            datamatrix[16][scounter]= sstdvec[3]

        scounter = scounter + 1 
        if(len(bar_data_raw)):
            datamatrix[7][0] = np.mean(bar_data_raw)
            datamatrix[8][0] =  np.std(bar_data_raw)
            for pos in range(0,4):
                datamatrix[17][pos] = np.mean(crosslinker_n_fil_bar[pos])
                datamatrix[18][pos] =  np.std(crosslinker_n_fil_bar[pos])
        
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker∂∂
    pd.DataFrame(datamatrix, index=['Time', 'Ncross_mean','Ncross_std','Valency_mean','Valency_std',
                                    'Nfil_per_cross_mean','Nfil_per_cross_std',
                                    'Bar_Nfil_per_cross_mean','Bar_Nfil_per_cross_std',
                                    'Mean_Nfpc_1_fil','Std_Npc_1_fil',
                                    'Mean_Nfpc_2_fil','Std_Nfpc_2_fil',
                                    'Mean_Nfpc_3_fil','Std_Nfpc_3_fil',
                                    'Mean_Nfpc_4_fil','Std_Nfpc_4_fil',
                                    'Bar_Mean_Nfpc','Bar_Std_Nfpc'
                                    ]).to_csv(outputfile)
    
def getvaspprops_MOVAVG(rset, deltasnap, N, Rval,outputfile):
    MEANID = 0
    STDID = 1
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    bar_min_snap = int(0.95*Nsnaps)
    bar_data_raw = np.array([])
    collectstatus = False
    datawritestatus = False
    datamatrix = np.zeros((9,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    # Nvasp = np.zeros((2,Ndatapoints))
    # Nvalency = np.zeros((2,Ndatapoints))
    # Nfil_per_solid = np.zeros((2,Ndatapoints))
    scounter = 0
    binsize = int(deltasnap/2)
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        nvasp_tmp = np.array([])
        nvalency_tmp = np.array([])
        Nfil_per_solid_tmp = np.array([])
        minsnap = max(SREF-binsize,0)
        maxsnap = min(SREF+binsize, Nsnaps)
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            for MSREF in range(minsnap,maxsnap):
                # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
                if(MSREF<=nsnaplocal):
                    # The two variables below together represent information of each bond in the network.
                    crossboundfilid = rset[runidx].snap[MSREF].crossboundfilid
                    crossboundblobid = rset[runidx].snap[MSREF].crossboundblobid
                    unique, counts = np.unique(crossboundblobid, return_counts=True)
                    nvasp_tmp = np.append(nvasp_tmp,len(unique))
                    nvalency_tmp = np.append(nvalency_tmp,np.array(counts))
                    # Get the filament IDs for each solid and calculate Nhands/Nfil
                    # This helps you understand the number of hands in each solid that are bound to the same filament.
                    # The distribution of Nhands/Nfil - closer to 1 suggests that each bound hand is bound to a different filament
                    # >1 suggests that multiple hands of the solid are bound to the same fil ID
                    # <1 is not possible.
                    for uiter, ublobid in enumerate(unique):
                        #Find the locs using argwhere
                        locs = np.argwhere(crossboundblobid==ublobid)
                        filid_solid = crossboundfilid[locs]
                        unique_fid, counts_fil = np.unique(filid_solid, return_counts=True)
                        lval = len(unique_fid)
                        Nfil_per_solid_tmp = np.append(Nfil_per_solid_tmp,lval)
        if(len(nvasp_tmp)):
            datamatrix[1][scounter]= np.mean(nvasp_tmp)
            datamatrix[2][scounter]= np.std(nvasp_tmp)
        if(len(nvalency_tmp)):
            datamatrix[3][scounter]= np.mean(nvalency_tmp)
            datamatrix[4][scounter]= np.std(nvalency_tmp)
        if(len(Nfil_per_solid_tmp)):
            datamatrix[5][scounter]= np.mean(Nfil_per_solid_tmp)
            datamatrix[6][scounter]= np.std(Nfil_per_solid_tmp)
        scounter = scounter + 1 
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(datamatrix, index=['Time', 'Ncross_mean','Ncross_std','Valency_mean','Valency_std',
                                    'Nfil_per_cross_mean','Nfil_per_cross_std',
                                    'Bar_Nfil_per_cross_mean','Bar_Nfil_per_cross_std']).to_csv(outputfile)
    


def getfillengthprops(rset, deltasnap, N, Rval, outputfile,outputfile2):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    print(Nsnaps,flush=True);
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    datamatrix = np.zeros((7,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    scounter = 0
    movsnapcounter =0
    Lmatrix = []
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        Nfilsnap = []
        L = []
        Lsum = []

        for runidx in range(0,len(rset)):
            Lsum_r = 0
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<nsnaplocal):
                # print('Shape of Nfilsnap='+str(np.shape(Nfilsnap)))
                # #print(rset[runidx].snap[SREF].filcoord)
                # print('Shape of filcoord='+str(len(rset[runidx].snap[SREF].filcoord)))
                # temp = rset[runidx].snap[SREF].filcoord
                # for f, fc in enumerate(temp):
                #     print(np.shape(fc));
                filcoord = (rset[runidx].snap[SREF].filcoord)
                Nfilsnap.append(len(filcoord))
                for f, fc in enumerate(filcoord):
                    Lfil_r = getfillength(fc)
                    Lsum_r = Lsum_r + Lfil_r
                    L.append(Lfil_r)
            Lsum.append(np.sum(Lsum_r))
        Lmatrix.append(L)
        datamatrix[1][scounter]= np.mean(Nfilsnap)
        datamatrix[2][scounter]= np.std(Nfilsnap)
        datamatrix[3][scounter]= np.mean(L)
        datamatrix[4][scounter]= np.std(L)
        datamatrix[5][scounter]= np.mean(Lsum)
        datamatrix[6][scounter]= np.std(Lsum)
        scounter = scounter + 1 
        
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker∂∂
    pd.DataFrame(datamatrix, index=['Time', 'Nfil_mean','Nfil_std','Lfil_mean','Lfil_std',
                                    'Lsum_mean','Lsum_std']).to_csv(outputfile)
    # Calculate moving average
    MOVAVGSCOUNTER = 0
    MOVAVG_NSNAPS = 5
    scounter = 0
    # This matrix will hold the PDF fillength
    #6 bins - 1micron, 12bins = 0.5, 24 - 0.25, 48,0.125
    fillen_bin_edges = np.linspace(0,6.5,100)
    #print(fillen_bin_edges)
    datamatrix2 = np.zeros((Ndatapoints,len(fillen_bin_edges)-1))
    for SREF in range(0,Nsnaps,deltasnap):
        minsnap = np.max([0,SREF-MOVAVG_NSNAPS])
        maxsnap = np.min([Nsnaps,SREF+MOVAVG_NSNAPS])
        Ltemp = []
        for MSREF in range(minsnap,maxsnap):
            Ltemp = Ltemp + Lmatrix[MSREF]
        pdfvec, bins = np.histogram(Ltemp,bins=fillen_bin_edges,density=True)
        datamatrix2[scounter][:] = pdfvec
        scounter = scounter + 1
    pd.DataFrame(datamatrix2, index=np.arange(0,Nsnaps,deltasnap)).to_csv(outputfile2)
    

def getradialdensitysnap(filcoord, Radialbinvec, histcountsvec):
    for f, fc in enumerate(filcoord):
        # Interpolate each monomer in the filament
        interpcoord = ic.interpolateallmonomers(fc)
        #interpcoord = fc
        Nbeads = (np.shape(interpcoord))[0]
        # Distance of each point from the center
        normvec = np.linalg.norm(interpcoord,axis=1)

        hist, bin_edges  = np.histogram(normvec,bins=Radialbinvec,density=False)
        histcountsvec = histcountsvec + hist
    return histcountsvec

def getradialdensitysnapVASP(filcoord, Radialbinvec, histcountsvec):
    print(filcoord is None)
    Nbeads = (np.shape(filcoord))[0]
    # Distance of each point from the center
    print(np.shape(filcoord))
    print(filcoord)
    normvec = np.linalg.norm(filcoord,axis=1)

    hist, bin_edges  = np.histogram(normvec,bins=Radialbinvec,density=False)
    histcountsvec = histcountsvec + hist
    return histcountsvec

def getradialdensity(rset, deltasnap, N, Rval, outputfile):
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    # Row represents data 
    datamatrix = np.zeros((5,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    densitydf = pd.DataFrame()
    Radialbinvec = np.linspace(0,1,41)
    Radialbincenter = Radialbinvec[0:len(Radialbinvec)-1]+(Radialbinvec[1]-Radialbinvec[0])
    densitydf["Radius"] = Radialbincenter
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        densitysnap = np.zeros(len(Radialbinvec)-1)
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        Nfilsnap = []
        L = []
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                filcoord = (rset[runidx].snap[SREF].filcoord)
                densitysnap = getradialdensitysnap(filcoord,Radialbinvec,densitysnap)
        densitydf[SREF] = densitysnap
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker∂∂
    densitydf.to_csv(outputfile)

def getradialdensityVASP(rset, deltasnap, N, Rval, outputfile):
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    # Row represents data 
    datamatrix = np.zeros((5,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    densitydf = pd.DataFrame()
    Radialbinvec = np.linspace(0,1,41)
    Radialbincenter = Radialbinvec[0:len(Radialbinvec)-1]+(Radialbinvec[1]-Radialbinvec[0])
    densitydf["Radius"] = Radialbincenter
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        densitysnap = np.zeros(len(Radialbinvec)-1)
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        Nfilsnap = []
        L = []
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                solidcoord = np.array(rset[runidx].snap[SREF].solidcoord)
                densitysnap = getradialdensitysnapVASP(solidcoord,Radialbinvec,densitysnap)
        densitydf[SREF] = densitysnap
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker∂∂
    densitydf.to_csv(outputfile)

def geteigensnap(filcoord):
    snapcoordmat=np.array([[]])
    for f, fc in enumerate(filcoord):
        Nbeads = (np.shape(fc))[0]
        filcoordmat = np.zeros((Nbeads,3))
        for i in range(0,Nbeads):
            filcoordmat[i,:] = np.array(fc[i,:])
        if f==0:
            snapcoordmat = filcoordmat.copy()
        else:
            axisval = 0
            snapcoordmat = np.concatenate((snapcoordmat,filcoordmat),axis=axisval)
    
    COM = np.mean(snapcoordmat,0)
    # Calculate PC1
    snapcoordmat = (snapcoordmat - COM)
    crosscorr  = np.matmul(np.transpose(snapcoordmat),snapcoordmat)/len(snapcoordmat)
    eval, evec =np.linalg.eig(crosscorr)
    eval = np.flip(np.sort(eval))
    eval = 2*np.sqrt(eval)
    return eval

def geteigentrajectory(rset, deltasnap, N, Rval, outputfile):
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    # Row represents data 
    datamatrix = np.zeros((4,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        Nfilsnap = []
        L = []
        PC1 = [];PC2=[];PC3=[];
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                filcoord = (rset[runidx].snap[SREF].filcoord)
                PCeigen = geteigensnap (filcoord)
                PC1.append(PCeigen[0])
                PC2.append(PCeigen[1])
                PC3.append(PCeigen[2])
        datamatrix[1][scounter] = np.mean(PC1)
        datamatrix[2][scounter] = np.mean(PC2)
        datamatrix[3][scounter] = np.mean(PC3)
        scounter = scounter + 1
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker
    pd.DataFrame(datamatrix,['Time','eig1','eig2','eig3']).to_csv(outputfile)


def analyzetrajectory(fpathvar,N, dirname, outfilename):
    f=4
    print(dirname,flush=True)
    meshvar=ic.icosphere(f,1.0);
    vertices = meshvar[0]
    triangles = meshvar[1]
    Ntriangles = np.shape(triangles)[0]
    unitnormal = ic.generateunitnormals(meshvar,1.0)
    tarea = ic.gettrianglearea(meshvar)
    Nlists = ic.generateNeighborList(meshvar)
    print('mesh created', flush=True)
    #Generate
    deltasnap = 1
    Nruns = 16
    Nreps = 5
    Ntotal = Nruns*Nreps
    Nsnaps = 601
    flist=np.zeros((1,Nsnaps))
    counter = 0
    dirlist = []
    dirlist.append(dirname)
    r = readcytosimtraj(fpathvar+dirname+'/')
    getfracoccupiedtimeseries(N, dirname, meshvar, Nlists, unitnormal, Ntriangles, r, deltasnap, outfilename)
    
def frontend_FO(N, Nreps):
    if N==1:
        foldername = 'Fgrow_10_3_N_1e3_5reps'
    else:
        foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
    fpathvar = '/scratch/achndrasekaran/cytosimfiles/'+foldername+'/'    
    for ridx in range(0,nr):
        for cidx in range(0,nc):
            Rval = nc*ridx+cidx
            for repid in range(0,Nreps):
                dirname = 'R_'+str(Rval)+'_r_'+str(repid)
                outfilename = foldername+'_'+dirname
                analyzetrajectory(fpathvar,N, dirname, outfilename)
    print("The End....", flush=True)

def frontend_FO_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    if(len(args)==0):
        outputfile = 'FO_set_'+str(N)+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/FO_'+foldername+'_'+dirname+'.csv'
    print(foldername, flush=True)
    analyzetrajectory(fpathvar,N, dirname, outputfile)
    print("The End....", flush=True)

def frontend_cross_set_Rval(N, Rval, Nreps, *args):
    if(len(args)==0):
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
        outputfile = 'Crosslink_set_'+str(N)+'_R_'+str(Rval)+'.csv'
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Crosslink_set_'+str(N)+'_R_'+str(Rval)++'.csv'
    print(foldername, flush=True)
    rset = []
    for repid in range(0,Nreps):
        dirname = 'R_'+str(Rval)+'_r_'+str(repid)
        print('Reading trajectory '+dirname,flush=True)
        r = readcytosimtraj(fpathvar+dirname+'/')
        print('Trajectory read..',flush=True)
        rset.append(r[0])
        print('Read repid '+str(repid),flush=True)
    getvaspprops(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_cross_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    if(len(args)==0):   
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'
        outputfile = 'Crosslink_set_'+str(N)+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    else:
        print(args, flush=True)
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getvaspprops(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_cross_set_Rval_MOVAVG(N, Rval, Nreps):
    if N==1:
        foldername = 'Fgrow_10_3_N_1e3_5reps'
    else:
        foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
    print(foldername, flush=True)
    fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'    
    rset = []
    for repid in range(0,Nreps):
        dirname = 'R_'+str(Rval)+'_r_'+str(repid)
        print('Reading trajectory '+dirname,flush=True)
        r = readcytosimtraj(fpathvar+dirname+'/')
        print('Trajectory read..',flush=True)
        rset.append(r[0])
        print('Read repid '+str(repid),flush=True)
    outputfile = 'Crosslink_MA_set_'+str(N)+'_R_'+str(Rval)+'.csv'
    getvaspprops_MOVAVG(rset, 5, N, Rval, outputfile)
    print("The End....", flush=True)

def frontend_mean_std_FO(N, Rval, Nreps):
    if N==1:
        foldername = 'Fgrow_10_3_N_1e3_5reps'
    else:
        foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
    fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/' 

def frontend_filprops_Rval(N, Rval, Nreps, *args):
    if(len(args)==0):
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
        outputfile = 'Fillength_set_'+str(N)+'_R_'+str(Rval)+'.csv'
        outputfile2 = 'Fillength_set_dist'+str(N)+'_R_'+str(Rval)+'.csv'
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Fillength_'+foldername+'R_'+str(Rval)+'.csv'
        outputfile2 = args[1]+'outputfiles/Fillength_dist_'+foldername+'R_'+str(Rval)+'.csv'
    print(foldername, flush=True)
    rset = []
    for repid in range(0,Nreps):
        dirname = 'R_'+str(Rval)+'_r_'+str(repid)
        print('Reading trajectory '+dirname,flush=True)
        r = readcytosimtraj(fpathvar+dirname+'/','filonly')
        print('Trajectory read..',flush=True)
        rset.append(r[0])
        print('Read repid '+str(repid),flush=True)
    getfillengthprops(rset, 1, N, Rval,outputfile, outputfile2)
    print("The End....", flush=True)

def frontend_density_Rval_repid(N, Rval, repid, *args):
    if(len(args)==0):
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
        outputfile = 'Fillength_set_'+str(N)+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Density_'+foldername+'R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)
    rset = []
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/','filonly')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getradialdensity(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_densityVASP_Rval_repid(N, Rval, repid, *args):
    if(len(args)==0):
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
        outputfile = 'Fillength_set_'+str(N)+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Density_'+foldername+'R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)
    rset = []
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/','filonly')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getradialdensityVASP(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_eigval_Rval_repid(N, Rval, repid, *args):
    if(len(args)==0):
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
        outputfile = 'Eigval_set_'+str(N)+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Eigval'+foldername+'R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)
    rset = []
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/','filonly')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    geteigentrajectory(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_density_Rval(N, Rval, Nreps, *args):
    if(len(args)==0):
        if N==1:
            foldername = 'Fgrow_10_3_N_1e3_5reps'
        else:
            foldername = 'Fgrow_10_3_N_1e3_5reps_set'+str(N)
        fpathvar = '/scratch/achandrasekaran/cytosimfiles/'+foldername+'/'   
        outputfile = 'Fillength_set_'+str(N)+'_R_'+str(Rval)+'.csv'
    else:
        foldername = args[0]
        fpathvar = args[1]+'/'+foldername+'/'  
        outputfile = args[1]+'outputfiles/Density_'+foldername+'R_'+str(Rval)+'.csv'
    print(foldername, flush=True)
    rset = []
    for repid in range(0,Nreps):
        dirname = 'R_'+str(Rval)+'_r_'+str(repid)
        print('Reading trajectory '+dirname,flush=True)
        r = readcytosimtraj(fpathvar+dirname+'/','filonly')
        print('Trajectory read..',flush=True)
        rset.append(r[0])
        print('Read repid '+str(repid),flush=True)
    getradialdensity(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)