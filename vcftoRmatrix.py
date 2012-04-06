import simplejson

class genotypes:
    def __init__(self, name):
        self.name = name
        self.snpfile = open('./' + name + 'SnpPos')
        self.linesfile = open('./' + name + 'Lines')
        self.genofile = open('./' + name + 'Geno')
        self.output = open('./'+name+'output', 'w')

    def __selectSNPs__(self, chosenSNPs):
        snps = simplejson.load(self.snpfile)
        self.csindex = map(lambda x: snps.index(x), filter(lambda y: y in snps, chosenSNPs))



def combineflatgenos(names, chosenSNPs):
    files  = map(lambda x: genotypes(x), names)
    for f in files:
        f.__selectSNPs__(chosenSNPs)
        for c in f.csindex:
            f.genofile.seek((x-1)*59*2)
            r = f.genofile.read(59*2)
            print(r)
            f.output.write(r)
            f.genofile.seek(0)


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
                

# turn vcf file into a hash (python dictionary)
def hashvcf(vcffile):
    file = open(vcffile)
    hash = {}
    lines = file.readlines(10000000)
    while(lines!= []):
        for l in lines:
            if not l.startswith('#'):
                tokens = l.strip('\n').split('\t')
                key = 'chr'+tokens[0]+'pos'+tokens[1]
                hash[key] = []
                for t in tokens[9:]:
                    if t =='.':
                        matrixentry = 0
                    else:
                        matrixentry = int(t[0])+int(t[2])
                        hash[key].append(matrixentry)
        lines = file.readlines(10000000)
    return hash

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

def convertIlluminaSNPstoBED():
    file = open('./Human1M-SNPlist.txt')
    lines = file.readlines()
    outputfile = open('./IlluminaHG19BED','w')
    lines = lines[1:]
    for l in lines:
        t = l.split('\t')
        outputfile.write('chr'+t[2] + '\t'+t[3].strip('\n').strip('\r')+ '\t' + str(int(t[3]) + 1) + '\n')


