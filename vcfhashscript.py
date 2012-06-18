import simplejson
import vcftoRmatrix 
import threading
from threading import Thread

#vcftoRmatrix.flatfilevcf('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../genotypes/CEUlowcov')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../genotypes/YRIlowcov')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../genotypes/CHBJPTlowcov')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../genotypes/YRItrio')
#vcftoRmatrix.flatfilevcf('../1000GenomesData/CEU.trio.2010_09.genotypes.vcf', '../genotypes/CEUtrio')


file = open('./omnisnpsHASH')
omnisnps = simplejson.load(file)
file.close()
snps = omnisnps.keys()

i = 0    
size = 1000
while(i*size < len(snps)):
    if (i+1) * size < len(snps):
        smpl = snps[i*size : (i+1)*size]
    else:
        smpl = snps[i*size:]
    Thread(target = vcftoRmatrix.combineflatgenos(['../genotypes/CEUlowcov','../genotypes/CHBJPTlowcov', '../genotypes/YRIlowcov', '../genotypes/CEUtrio', '../genotypes/YRItrio'], smpl,i)).start()
    i=i+1
    
#files = vcftoRmatrix.combineflatgenos(['../genotypes/CEUlowcov','../genotypes/CHBJPTlowcov', '../genotypes/YRIlowcov', '../genotypes/CEUtrio', '../genotypes/YRItrio'], snps)
#files[0]



