import os, sys

spp = ['Vc', 'Vm', 'Vr', 'Vv']
filelist = os.listdir('.')
#print filelist
buscos = {}
for sp in spp:
    buscos[sp] = []
for filename in filelist:
    for thing in spp:
        if thing in filename and 'completebuscos.fasta' in filename:
            sp = thing
    data = open(filename, 'r')
    for line in data:
        line = line.strip()
        if line.startswith('>'):
            ID = line[1:].split('|')[0]
            buscos[sp].append(ID)
    data.close()
qbuscos = []
allbuscos = []
for sp in spp:
    allbuscos = allbuscos + list(set(buscos[sp]))
for thing in set(allbuscos):
    if allbuscos.count(thing) == 4:
        qbuscos.append(thing)
print len(qbuscos)

print qbuscos[:5]

buscoseqs = {}
for sp in spp:
    fetched = 0
    buscoseqs[sp] = ['ok']*len(qbuscos)
    for filename in filelist:
        if sp in filename and 'completebuscos.fasta' in filename:
            data = open(filename, 'r')
            for line in data:
                line = line.rstrip()
                if line.startswith('>'):
                    ID = line[1:].split('|')[0]
                    if ID in qbuscos and buscoseqs[sp][qbuscos.index(ID)] == 'ok':
                        buscoseqs[sp][qbuscos.index(ID)] = data.next().strip()
                        fetched += 1
            data.close()
            print filename, fetched


dataset = {}
for sp in spp:
    dataset[sp] = ''

alndir = '/home/clone/NEWVEXILLUM/aligned_buscos/'
for thing in qbuscos:
    writefile = open(alndir+'/'+thing+'.fas', 'w')
    for key in buscoseqs:
        print >>writefile, '>'+key+'\n'+buscoseqs[key][qbuscos.index(thing)]
    writefile.close()
    cmd = 'mafft --thread 12 ' + alndir+'/'+thing+'.fas' + ' > ' + alndir+'/'+thing+'_aligned.fas'
    os.system(cmd)
    cmd2 = 'python '+alndir+'makesomethingNotInterleaved.py '+alndir+thing+'_aligned.fas '+alndir+thing+'_alignedNI.fas'
    os.system(cmd2)
    aligned = open(alndir+'/'+thing+'_alignedNI.fas', 'r')
    for line in aligned:
        line = line.rstrip()
        if line.startswith('>'):
            dataset[line[1:]] = dataset[line[1:]]+aligned.next().strip()
    aligned.close()

variable = 0
thelen = len(dataset.values()[0])
print thelen

for i in range(thelen):
    AAlist = [dataset[key][i] for key in dataset if dataset[key][i] != '-']
    if len(set(AAlist)) > 1:
        variable += 1
print variable
finalwrite = open('Vex_buscofasta.fas', 'w')
for key in dataset:
    print >>finalwrite, '>'+key+'\n'+dataset[key]
finalwrite.close()

iqtree = 'iqtree -nt 16 -st AA -s Vex_buscofasta.fas -m Dayhoff -bb 1000 -safe'
os.system(iqtree)



   
