import simplejson
import vcftoRmatrix 
import threading
from threading import Thread

#vcftoRmatrix.flatfilevcf('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../genotypes/CEUlowcov')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../genotypes/YRIlowcov')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../genotypes/CHBJPTlowcov')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../genotypes/YRItrio')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/CEU.trio.2010_09.genotypes.vcf', '../genotypes/CEUtrio')


#file = open('./omnisnpsHASH')
#omnisnps = simplejson.load(file)
#file.close()
#snps = omnisnps.keys()


names = ['../genotypes/CEUlowcov','../genotypes/CHBJPTlowcov', '../genotypes/YRIlowcov', '../genotypes/CEUtrio', '../genotypes/YRItrio']

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
    newgenofile = open('chosenlinesSfile','w')
    genolines = file.readlines()
    print inds
    for l in genolines:
	l = l.strip('\n')
	newl = ','.join(map(lambda x: l[x], inds))
	newgenofile.write(newl + '\n') 

test = ["NA06985", "NA06986", "NA06994", "NA07000", "NA07037", "NA07051", "NA07346"]
getlines(names, 'mergedoutput', test)    
	
    



