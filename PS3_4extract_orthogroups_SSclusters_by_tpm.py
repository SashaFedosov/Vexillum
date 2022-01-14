import sys, os
maxtpmvalue = 200
def maxTpms(spmtpm):
    recs = spmtpm.split('|')
    tpms = []
    for rec in recs:
        tpm = float(rec.split('-')[1])
        tpms.append(tpm)
    return max(tpms)


SSdir = '/home/clone/NEWVEXILLUM/CDHIT'
ORTHOdir = '/home/clone/NEWVEXILLUM/ORTHOFINDER/OrthoFinder/Results_Sep23/Orthogroup_Sequences'
homedir = '/home/clone/NEWVEXILLUM/ORTHOFINDER'
qOGs = '/home/clone/NEWVEXILLUM/ORTHOFINDER/HiExpOGs/'
#os.mkdir(qOGs)
clusters06 = []
clusters65 = []
OG = []
allOG = {}
datafile = open(homedir+'/VexCMRV.alldata', 'r')
datafile.readline()
for line in datafile:
    line = line.strip()
    data = line.split('\t')
    if len(data) != 18:
        print len(data)
    else:
        tpm = data[6]
        if maxTpms(tpm) >= maxtpmvalue:
            clusters06.append(data[4])
            clusters65.append(data[5])
            OG.append(data[1])
            if data[1] in allOG.keys():
                allOG[data[1]]+=1
            else:
                allOG[data[1]]=1
datafile.close()

print len(set(clusters06)), len(set(clusters65)), len(set(OG))

ORFstoprint = 0
newdatafile = open(homedir+'/VexCMRV_NEWclusters_tpm'+str(maxtpmvalue)+'.alldata', 'w')
datafile1 = open(homedir+'/VexCMRV.alldata', 'r')
headers = datafile1.readline()
print >>newdatafile, headers
for line in datafile1:
    line=line.strip()
    data = line.split('\t')
    if len(data) != 18:
        print len(data)
    else:
        if data[1] in OG or data[4] in clusters06 or data[5] in clusters65:
            ORFstoprint += 1
            print >>newdatafile, line
            cmd = 'cp '+ORTHOdir+'/'+data[1]+'.fa '+qOGs
            os.system(cmd)
datafile.close()
newdatafile.close()     
print ORFstoprint