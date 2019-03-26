import os, glob, sys
import numpy as np
import pandas as pd
import csv

class rowObject:
    def __init__(self, row):
        self.name=row[0]
        self.segment=row[1]
        self.ntpos=int(row[2])
        self.binocheck=row[3]
        self.major=row[4]
        self.majorfreq=float(row[5])
        self.minor=row[6]
        self.minorfreq=float(row[7])
        self.afreq=float(row[8])
        self.cfreq=float(row[9])
        self.gfreq=float(row[10])
        self.tfreq=float(row[11])
        self.totalcount=int(float(row[12]))
        self.ref_nt=row[13]
        self.conmajorcheck=row[14]
        self.aa_pos=row[15]
        self.major_codon=row[16]
        self.major_aa=row[17]
        self.minor_codon=row[18]
        self.minor_aa=row[19]
        self.nonsyn_syn=row[20]

def l1_norm(somestrain, CUTOFF = 0.03, COVERCUTOFF = 200):

    path = 'FILES/fullvarlist/'
    count = 0
    STRAIN = somestrain.upper()
    for SEGMENT in ['HA']:#refdict:
        SEGMENT = STRAIN +'-'+ SEGMENT
        keylist = []
        VARDICT = {}
        for infile in glob.glob( os.path.join(path, '*'+STRAIN+'*'+SEGMENT+'.0.01.snplist.csv') ):
            rootname = infile.split('/')[-1].split('.')[0]
            df = pd.read_csv(infile)
            VARDICT[rootname] = df

            keylist.append(rootname)
        print keylist

        dismatrix = np.zeros((len(keylist),len(keylist)))

        for aidx,asamp in enumerate(keylist): #column
            for bidx,bsamp in enumerate(keylist): #row
                
                print asamp, bsamp
                if aidx == bidx:
                    dismatrix[aidx,bidx] = 0
                elif aidx > bidx: #Because matrix is symmetrical we cut work in half
                    dismatrix[aidx,bidx] = dismatrix[bidx,aidx]
                else:
                    adf = VARDICT[asamp] 
                    bdf = VARDICT[bsamp]
                    alist = []
                    blist = []
                    for idx,arow in adf.iterrows():
                        if arow['totalcount'] == 0:
                            alist.append([0,0,0,0.0])
                        else:
                            avector = [arow['A']/float(arow['totalcount']),arow['C']/float(arow['totalcount']),arow['G']/float(arow['totalcount']),arow['T']/float(arow['totalcount'])]
                            alist.append(avector)
                    for idx,arow in bdf.iterrows():
                        if arow['totalcount'] == 0:
                            blist.append([0,0,0,0.0])
                        else:
                            avector = [arow['A']/float(arow['totalcount']),arow['C']/float(arow['totalcount']),arow['G']/float(arow['totalcount']),arow['T']/float(arow['totalcount'])]
                            blist.append(avector)
                    val = 0
                    for a,b in zip(alist,blist):
                        a = np.array(a)
                        b= np.array(b)
                        val = val + np.linalg.norm((a-b),ord=1)
                    dismatrix[aidx,bidx] = val

        thefile = open(somestrain+'.namelist.csv','w')
        for key in keylist:
            print>>thefile,key
        with open(somestrain+'.'+SEGMENT+'.dissim_all.key.csv', 'w') as f:
            csvWriter = csv.writer(f)
            csvWriter.writerow(keylist)
            for idx,row in enumerate(np.atleast_2d(dismatrix.T)):
                csvWriter.writerow(row)

l1_norm('H1N1')
l1_norm('H3N2')



