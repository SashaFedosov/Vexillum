import sys, os

spp = ['Vc', 'Vm', 'Vr', 'Vv']

ORFdir = '/home/clone/NEWVEXILLUM/ORFs/'
depthdir = '/home/clone/NEWVEXILLUM/clustering/'

ORFlist = []
ORFfile = open(depthdir + 'all_decont_ORFs.txt', 'r')
ORFfile.readline()
for line in ORFfile:
    line = line.rstrip()
    ORFlist.append(line.split('\t')[0])
ORFfile.close()
print len(ORFlist)

alldepths = open(depthdir+'all_decontORFs_depthdata.txt', 'w')

misassemblies = {}
depthdata = {}
ORFdata = {}
for sp in spp:
    spms = []
    sporfs = [thing for thing in ORFlist if sp in thing]
    datafile = open(ORFdir+sp+'.alldata', 'r')
    for line in datafile:
        line = line.rstrip()
        data = line.split('\t')
        if data[0] in sporfs:
            orfID = data[2]
            if 'Contig' in orfID:
                contID = orfID.split('_')[1].split(':')[0]
            else:
                contID = 'TRINITY' + orfID.split('TRINITY')[1].split(':')[0]
            start = int(orfID.split(':')[1])
            end = int(orfID.split(':')[2])
            spm = data[1].replace('39', '35')
            if not spm in spms:
                spms.append(spm)
            TPM = float(data[3])
            if data[0] in ORFdata.keys() and ORFdata[data[0]][-1] >= TPM:
                continue
            else:
                ORFdata[data[0]] = [spm, contID, min([start, end]), max([start, end]), TPM]
                depthdata[data[0]] = []
    datafile.close()
    print sp, len(sporfs)
    
    for spm in spms:
        spmORFs = []
        spmconts = []
        depths = {}
        for thing in sporfs:
            if ORFdata[thing][0] == spm:
                spmORFs.append(thing)
                spmconts.append(ORFdata[thing][1])
        print spm, len(spmORFs), len(spmconts)
        print spmORFs[0], spmconts[0]
        
        depthfile = open(depthdir + spm + '_depth.out', 'r')
        firstline = depthfile.readline()
        CONT = firstline.split('\t')[0]
        depthfile.close()
        print CONT
        
        depthfile = open(depthdir + spm + '_depth.out', 'r')
        good = False
        for line in depthfile:
            line = line.strip()
            Cont = line.split('\t')[0]
            site = int(line.split('\t')[1])
            depth = int(line.split('\t')[2])
            if Cont == CONT:
                if good == True:
                    depths[Cont].append([site, depth])
            else:
                CONT = Cont
                if Cont in set(spmconts):
                    depths[Cont] = [[site, depth]]
                    good = True
                else:
                    good = False
        depthfile.close()
        print len(depths)
        
        if spm == '8Vcsg':
            print 'TRINITY_DN399_c0_g1_i3' in depths.keys()
        
        for thing in spmORFs:
            if ORFdata[thing][1] not in depths.keys():
                print thing, ORFdata[thing][1], 'MISSING'
            else:
                depthrecords = depths[ORFdata[thing][1]]
                OK = True
                for item in depthrecords:
                    if item[0] >= ORFdata[thing][2] and item[0] <= ORFdata[thing][3]:
                        depthdata[thing].append(item[1])
                for i in range(len(depthdata[thing])-1):
                    if depthdata[thing][i]*2 < depthdata[thing][i+1] and depthdata[thing][i] >= 10:
                        misassemblies[thing] = depthdata[thing]
                        print thing, depthdata[thing][i], depthdata[thing][i+1]
                        OK = False
                    elif depthdata[thing][i] > depthdata[thing][i+1]*2 and depthdata[thing][i+1] >= 10:
                        misassemblies[thing] = depthdata[thing]
                        print thing, depthdata[thing][i], depthdata[thing][i+1]
                        OK = False
                    else:
                        continue
                if OK == True:
                    print >>alldepths, thing+'\t'+'OK'+'\t'+str(depthdata[thing])
                else:
                    print >>alldepths, thing+'\t'+'MISASSEMBLED'+'\t'+str(depthdata[thing])
alldepths.close()
        
        
    
    
    
