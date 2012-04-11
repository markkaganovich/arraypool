import simplejson

class genotypes:
    def __init__(self, name):
        self.name = name
        self.snpfile = open('./' + name + 'SnpPos')
        self.linesfile = open('./' + name + 'Lines')
        self.genofile = open('./' + name + 'Geno')
        self.tempoutput = open('./'+name+'tempoutput', 'w')
        self.__genolen__()

    def __selectSNPs__(self, chosenSNPs):
        self.snps = simplejson.load(self.snpfile)
        self.csindex = map(lambda x: self.snps.index(x), filter(lambda y: y in self.snps, chosenSNPs))
        self.havechosensnps = map(lambda x: self.snps[x], self.csindex)

    def __genolen__(self):
        l = self.genofile.readline()
        self.ln = len(l.split(','))
        self.genofile.seek(0)
        
def combineflatgenos(names, chosenSNPs):
    files  = map(lambda x: genotypes(x), names)
    for f in files:
        f.__selectSNPs__(chosenSNPs)

    unionchosensnps = []
    map(lambda x: unionchosensnps.extend(x.havechosensnps), files) 
    unionchosensnps = list(set(unionchosensnps))
   
    output = open('ceumergetest','w')
    for s in unionchosensnps:
        totalr = ''
        for f in files:
            if s in f.havechosensnps:
                c = f.snps.index(s)
                f.genofile.seek((c-1)*f.ln*2)
                r = f.genofile.read(f.ln*2)
                f.genofile.seek(0)
            else:
                r = '0,' * f.ln
            totalr = totalr +r.strip('\n')
        output.write(totalr.strip(',')+'\n')
                
    output.close()

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


