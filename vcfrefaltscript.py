import simplejson
import vcftoRmatrix

vcftoRmatrix.getrefalt('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../genotypes/19CEUlowcov')
vcftoRmatrix.getrefalt('../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../genotypes/19YRIlowcov')
vcftoRmatrix.getrefalt('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../genotypes/19CHBJPTlowcov')
vcftoRmatrix.getrefalt('../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../genotypes/19YRItrio')
vcftoRmatrix.getrefalt('../1000GenomesData/CEU.trio.2010_09.genotypes.vcf', '../genotypes/19CEUtrio')
