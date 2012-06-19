import simplejson

class genotypes:
    def __init__(self, name):
        self.name = name
        self.snpfile = open('./' + name + 'SnpPos')
        self.linesfile = open('./' + name + 'Lines')
        self.genofile = open('./' + name + 'Geno')
        self.__genolen__()
    
    def __selectSNPs__(self, chosenSNPs):
        self.chosenoutput = open('./' + self.name + 'chosenoutput' + str(mapnum), 'w')
        self.snps = simplejson.load(self.snpfile)
        self.csindex = []
        for x in chosenSNPs:
			if x in self.snps:
				self.csindex.append(self.snps.index(x))
			else:
				self.csindex.append('NA')
			#self.havechosensnps = map(lambda x: self.snps[x], self.csindex)
        del self.snps
    
    def __genolen__(self):
        l = self.genofile.readline()
        self.ln = len(l.split(','))
        self.genofile.seek(0)

def convertprecombine(files, chosenSNPs):
    for f in files:
        f.__selectSNPs__(chosenSNPs)
        writeprecombinefile(f, mapnum)
        
        #Thread
        i = 0
        size = 1000
        while(i*size < len(chosenSNPs)):
            if (i+1) * size < len(chosenSNPs)
                smpl = snps[i*size : (i+1)*size]
            else:
                smp= snps[i*size:]
            Thread(target = writeprecombinefile(f, i)).start()
            i = i+1


def writeprecombinefile(f, mapnum)
        f.chosenoutput = open('./' + self.name + 'chosenoutout' + str(mapnum), 'w')
    	for c in f.csindex:
            if c != 'NA':
                f.genofile.seek((c-1)*f.ln*2)
                r = f.genofile.read(f.ln*2)
                f.genofile.seek(0)
                f.chosenoutput.write(r)
            else:
                r = '0,' * f.ln
                f.chosenoutput.write(r.strip(',')+'\n')
        f.chosenoutput.close()

def combineflatgenos(names, chosenSNPs):
    files = map(lambda x: genotypes(x), names)
    convertprecombine(files, chosenSNPs)
    file = open('mergedoutput','w')
    for f in files:
        f.opened = open('./' + f.name + 'chosenoutput')
    for i in range(0,len(chosenSNPs)):
        l = ''
        for f in files:    
            l = l + f.opened.readline().strip('\n') + ','
        file.write(l.strip(',')+'\n')


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
                snppos.append('chr'+tokens[0]+'pos'+tokens[1])
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


