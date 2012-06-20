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

vcftoRmatrix.combineflatgenos(['../genotypes/CEUlowcov','../genotypes/CHBJPTlowcov', '../genotypes/YRIlowcov', '../genotypes/CEUtrio', '../genotypes/YRItrio'], snps)
    



