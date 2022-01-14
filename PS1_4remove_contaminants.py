# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 06:40:46 2021

@author: Sasha
"""
import numpy as np
import os, sys

contamination_threshold = 0.01
Mazur = ['8VcgL', '8Vcsg', '9VcgL', '9Vcsg', '13VvgL', '13Vvsg', '14VvgL', '14Vvsg']
Scoltech = ['33VmgL', '33Vmsg', '35VmgL', '39Vmsg', '44VrsggL']
allspms = Mazur + Scoltech

OrfIDDict = {}
AAseqs = []
ORFtpms = {}
datafile = open('/home/clone/NEWVEXILLUM/FILTER_final_ORFs/VexCMRV_NEWclusters_tpm200.alldata', 'r')
datafile.readline()
datafile.readline()
for line in datafile:
    Line = line.rstrip()
    data = line.split('\t')
    OrfIDDict[data[0]] = [data[6], data[8]]
    if not data[8] in AAseqs:
        AAseqs.append(data[8])
        ORFtpms[data[8]] = [0]*13
    sptpm = data[6].split('|')
    for thing in sptpm:
        spmrec = thing.split(':')[0]
        tpmrec = float(thing.split('-')[1])
        for i in range(len(allspms)):
            if allspms[i] == spmrec:
                ORFtpms[data[8]][i] += tpmrec
datafile.close()    
print len(OrfIDDict)
print len(AAseqs)
print allspms

print ORFtpms['MLNLTAQLAFFLACLALFSAVSATTCEDSSLKCAGWRLYCNKGNGVTDDIVTKARKMCPVTCGLCGSSECRDHLDECTSTLRDQCSKDVQVKKVFCRRTCGDCTPHGK*']

def ZeroDownContamins(maz):
    yeslist = ['yes']*len(maz)
    for k in range(len(maz)):
        if maz[k] == 0:
            yeslist[k] = 'no'
        else:      
            for m in range(len(maz)):
                if maz[k] < maz[m]*contamination_threshold:
                    yeslist[k] = 'no'
                    break
    for l in range(len(yeslist)):
        if yeslist[l] == 'no':
            maz[l] = 0
    return maz

def CompileTpm(spms, maz):
    if maz.count(0) == len(maz):
        tpmstring = 'X|'
    else:
        tpmstring = ''
        for i in range(len(maz)):
            if maz[i] != 0:
                tpmstring = tpmstring + spms[i] + ':TPM-'+str(maz[i])+'|'
    return tpmstring[:-1]
    

allnewtpms = {}
for seq in AAseqs:
    mazur = ORFtpms[seq][:8]
    scol = ORFtpms[seq][8:]
    newtpms = ZeroDownContamins(mazur) + ZeroDownContamins(scol)
    allnewtpms[seq] = newtpms
    if seq == 'MLNLTAQLAFFLACLALFSAVSATTCEDSSLKCAGWRLYCNKGNGVTDDIVTKARKMCPVTCGLCGSSECRDHLDECTSTLRDQCSKDVQVKKVFCRRTCGDCTPHGK*':
        print newtpms
Zeroed = 0

datafile = open('/home/clone/NEWVEXILLUM/FILTER_final_ORFs/VexCMRV_NEWclusters_tpm200.alldata', 'r')
headers = datafile.readline().strip()
datafile.readline()
newtpmfile = open('VexCMRV_NEWclusters_'+str(contamination_threshold*100)+'percent_decontaminated.txt', 'w')
print >>newtpmfile, headers.replace('Spms_and_tpms', 'Spms_and_tpms\tSpms_and_tpms_1percent decontaminated')
for line in datafile:
    line = line.rstrip()
    data = line.split('\t')
    ID = data[0]
    if ID.startswith('Vc'):
        newtpmrec = CompileTpm(allspms[:4], allnewtpms[data[8]][:4])
    elif ID.startswith('Vv'):
        newtpmrec = CompileTpm(allspms[4:8], allnewtpms[data[8]][4:8])
    elif ID.startswith('Vm'):
        newtpmrec = CompileTpm(allspms[8:12], allnewtpms[data[8]][8:12])
    else:
        newtpmrec = CompileTpm([allspms[12]], [allnewtpms[data[8]][12]])
    newline = line.replace(data[6], data[6]+'\t'+newtpmrec)
    print >>newtpmfile, newline
    if newtpmrec == 'X':
        Zeroed += 1
datafile.close()
newtpmfile.close()
print Zeroed





















