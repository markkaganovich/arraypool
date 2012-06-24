import simplejson
import threading
from threading import Thread

class genotypes:
    def __init__(self, name):
        self.name = name
        self.snpfile = open('./' + name + 'SnpPos')
        self.linesfile = open('./' + name + 'Lines')
        self.genofile = open('./' + name + 'Geno')
        self.__genolen__()
    
    def __genolen__(self):
        l = self.genofile.readline()
        self.ln = len(l.split(','))
        self.genofile.seek(0)

def calccsindex(filestruc, chosenSNPs):
    filestruc.csindex = []
    snps = simplejson.load(filestruc.snpfile)
    j = 0
    for x in chosenSNPs:
	print x
	for s in range(j, len(snps)):
	    if x == snps[s]:
		filestruc.csindex.append(s)
		j=s
		break
	    if s == len(snps)-1:
		print "break"
	        filestruc.csindex.append('NA')
	    
def convertprecombine(files, chosenSNPs):
    for f in files:
        print f.name
        calccsindex(f, chosenSNPs)
	writeprecombinefile(f)
	#Thread
        #i = 0
        #size = 100000
        #while(i*size < len(f.csindex)):
        #    if (i+1) * size < len(f.csindex):
        #        smpl = f.csindex[i*size : (i+1)*size]
        #    else:
        #        smpl = f.csindex[i*size:]
        #    Thread(target = writeprecombinefile(f, i, smpl)).start()
        #    i = i+size
	#    print i

def writeprecombinefile(f):
        chosenoutput = open('./' + f.name + 'outputprecombine', 'w')
    	for c in f.csindex:
            if c != 'NA':
                f.genofile.seek((c-1)*f.ln*2)
                r = f.genofile.read(f.ln*2)
                f.genofile.seek(0)
                chosenoutput.write(r)
            else:
                r = '0,' * f.ln
                chosenoutput.write(r.strip(',')+'\n')
        chosenoutput.close()

def combineflatgenos(names, chosenSNPs):
    files = map(lambda x: genotypes(x), names)
    convertprecombine(files, chosenSNPs)
    file = open('mergedoutput19','w')
    for f in files:
        f.opened = open('./' + f.name + 'outputprecombine')
    for i in range(0,len(chosenSNPs)):
        l = ''
        for f in files:    
            l = l + f.opened.readline().strip('\n') + ','
        file.write(l.strip(',')+'\n')

################################################################################################
#modifyied to the 19 version
def flatfilevcf(vcffile, outputname):
    file = open(vcffile)
    outputfile = open(outputname+'Geno', 'w')
    outputfileLines = open(outputname+'Lines', 'w')
    snppos = []
    outputfileSnpPos = open(outputname + 'SnpPos', 'w')
    lines = file.readlines(1000000)
    while(lines != []):
        for l in lines:
            if l.startswith('#CHROM'):
                genomenames = l.strip('\n').split('\t')
                simplejson.dump(genomenames[9:], outputfileLines)
                outputfileLines.close()
            if not l.startswith('#'):
                tokens = l.strip('\n').split('\t')
                f = filter(lambda x: 'GP' in x, tokens[7].split(';'))
		if f != []:
		    snppos.append('chr'+f[0].split('=')[1].split(':')[0]+'pos'+f[0].split('=')[1].split(':')[1])
                    m =''
                    for t in tokens[9:]:
                        m = m + str(int(t[0]) + int(t[2])) + ','
                    outputfile.write(m.strip(',')+'\n')
        lines = file.readlines(1000000)
    simplejson.dump(snppos, outputfileSnpPos)


def getvcfmatrix(filename, genomelist):
    file = open(filename)
    outputfile = open(filename+'matrixoutput', 'w')
    total = 0
    lines = file.readlines(1000000)
    while(lines!=[]):
        for line in lines:
            if not line.startswith('#'):
                total= total+1
                l = line.split('\t')
                for g in genomeind:
                    matrixentry = int(l[g][0]) + int(l[g][2])
                    outputfile.write(str(matrixentry) +'\t')
                outputfile.write('\n')
            elif line.startswith('#CHROM'):
                l = line.strip('\n').split('\t')
                genomeind = []
                for genome in genomelist:
                    genomeind.append(l.index('NA'+genome))
        lines = file.readlines(1000000)
        print total
    print total
    file.close()
    outputfile.close()


