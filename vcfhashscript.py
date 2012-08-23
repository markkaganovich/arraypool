import simplejson
import vcftoRmatrix 
import threading
from threading import Thread
'''
vcftoRmatrix.flatfilevcf('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../genotypes/19CEUlowcov')
print "1"
vcftoRmatrix.flatfilevcf('../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../genotypes/19YRIlowcov')
print "2"
vcftoRmatrix.flatfilevcf('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../genotypes/19CHBJPTlowcov')
vcftoRmatrix.flatfilevcf('../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../genotypes/19YRItrio')
vcftoRmatrix.flatfilevcf('../1000GenomesData/CEU.trio.2010_09.genotypes.vcf', '../genotypes/19CEUtrio')
'''
'''
file = open('./omni19snps')
snps = simplejson.load(file)
file.close()
'''
file = open('./correctedSNPs')
lines = file.readlines()
file.close()
snps = []
for l in lines:
    snps.append(l.strip('\n'))

names = ['../genotypes/19CEUlowcov','../genotypes/19CHBJPTlowcov', '../genotypes/19YRIlowcov', '../genotypes/19CEUtrio', '../genotypes/19YRItrio']

#vcftoRmatrix.combineflatgenos(names, snps)

def getlines(names, genotypefile, chosenlines):
    lines = []
    for n in names:
	LinesFile = open(n+'Lines')
	lines.extend(simplejson.load(LinesFile))
    inds = []
    for c in chosenlines:
	inds.append(lines.index(c))
    inds = map(lambda x: x *2, inds)
    file = open(genotypefile)
    #newgenofile = open('chosenlinesSfile19corrected','w')
    newgenofile = open('testlinesS','w')
    genolines = file.readlines()
    print inds
    for l in genolines:
	l = l.strip('\n')
	newl = ','.join(map(lambda x: l[x], inds))
	newgenofile.write(newl + '\n')
 
def getridofzero(genotypefile):
    lines = open(genotypefile).readlines()
    newoutputfile = open('notzero'+genotypefile, 'w')
    file = open(genotypefile +'nonzeroinds', 'w')
    for i in range(0, len(lines)):
	l = lines[i]
	x = map(lambda x: int(x), l.strip('\n').split(','))
        print x
	print l
	if sum(x) != 0:
	     print i
	     newoutputfile.write(l)
	     file.write(str(i) + '\n')

    file.close()
    
		
test = ["NA06985", "NA06986", "NA06994", "NA07000", "NA07037", "NA07051", "NA07346"]

file = open('testpool')
pool = simplejson.load(file)
getlines(names, 'mergedoutput19', pool)
#getridofzero('chosenlinesSfile19corrected')    	
getridofzero('testlinesS')    	
    



