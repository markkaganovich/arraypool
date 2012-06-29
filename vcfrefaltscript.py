import simplejson
import vcftoRmatrix

'''
vcftoRmatrix.getrefalt('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../genotypes/19CEUlowcov')
vcftoRmatrix.getrefalt('../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../genotypes/19YRIlowcov')
vcftoRmatrix.getrefalt('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../genotypes/19CHBJPTlowcov')
vcftoRmatrix.getrefalt('../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../genotypes/19YRItrio')
vcftoRmatrix.getrefalt('../1000GenomesData/CEU.trio.2010_09.genotypes.vcf', '../genotypes/19CEUtrio')


'''
names = ['../genotypes/19CEUlowcov','../genotypes/19CHBJPTlowcov', '../genotypes/19YRIlowcov', '../genotypes/19CEUtrio', '../genotypes/19YRItrio']

inputs = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']


for i in range(0,len(names)):
    name = names[i]
    input = inputs[i]
    vcftoRmatrix.flatfilevcf(input,name)

omniIDs = simplejson.load(open('./omni19IDs'))

omniIDsrefalt = {}

for name in names:
    snppos = simplejson.load(open(name+'RefAlt'))
    for id in omniIDs.keys():
        omniIDsrefalt[id] = snppos[omniIDs[id]]
    simplejson.dump(omniIDsrefalt, open(name+'RefAltIDs','w'))








