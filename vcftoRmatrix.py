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
	self.snpdict = open('./'+name+'SnpPosDic')    

    def __genolen__(self):
        l = self.genofile.readline()
        self.ln = len(l.split(','))
        self.genofile.seek(0)

def calccsindex(filestruc, chosenSNPs):
    filestruc.csindex = []
    snpdic = simplejson.load(filestruc.snpdict)
    inboth = set(snpdic.keys()) & set(chosenSNPs)
    for x in chosenSNPs:
	print x
	if x in inboth:
	    filestruc.csindex.append(snpdic[x])
	else:
	    filestruc.csindex.append('NA')
	    
def convertprecombine(files, chosenSNPs):
    for f in files:
        print f.name
        calccsindex(f, chosenSNPs)
	writeprecombinefile(f)

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

#add on
def getrefalt(vcffile, outputname):
    file = open(vcffile)
    outputfilerefalt = open(outputname+'RefAlt','w')
    lines = file.readlines(1000000)
    refalt=[]
    while(lines!=[]):
	for l in lines:
	    if not l.startswith('#'):
		tokens = l.strip('\n').split('\t')
		ref = tokens[3]
		alt = tokens[4]
		refalt.append({'ref':ref, 'alt':alt})
	lines = file.readlines(1000000)
    simplejson.dump(refalt, outputfilerefalt)





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


