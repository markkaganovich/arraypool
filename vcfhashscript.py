import simplejson
import vcftoRmatrix 

vcftoRmatrix.flatfilevcf('../1000Genomes/YRI.low_coverage.2010_09.genotypes.vcf', '../genotypes/YRIlowcov')
vcftoRmatrix.flatfilevcf('../1000Genomes/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../genotypes/CHBJPTlowcov')
vcftoRmatrix.flatfilevcf('../1000Genomes/YRI.trio.2010_09.genotypes.vcf', '../genotypes/YRItrio')
vcftoRmatrix.flatfilevcf('../1000Genomes/CEU.trio.2010_09.genotypes.vcf', '../genotypes/CEUtrio')


