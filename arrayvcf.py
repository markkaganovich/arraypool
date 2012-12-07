import vcf
import gl

poollines = ['NA19211', 'NA18943', 'NA19209', 'NA18526' ,'NA19161' , 'NA11920', 'NA11995' , 'NA18564' , 'NA18499' , 'NA12003']
names1KG = ['CEUlowcov','YRIlowcov','CHBJPTlowcov','YRItrio', 'CEUtrio']	
vcffiles = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf','../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', 
'../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']


def in1kg(poollines, names):
	lines = []
	for n in names:
		lines.extend(getlines(n))

def getlines(name):
	vfile = open(name, 'r')
	vcf_reader = vcf.Reader(vfile)
	lines = vcf_reader.samples
	return lines

#modified to the 19 version
def parse1KGvcf(vcffile, outputname):
	vfile = open(vcffile, 'r')
	vcf_reader = vcf.Reader(vfile)

	for record in vcf_reader:
		chrom = record.INFO['GP'].split(':')[0]
		pos = record.INFO['GP'].split(':')[1]	
		
	outputfile = open(outputname+'Geno', 'w')
	ref = {}
	alt = {}
	lines = file.readlines(1000000)
	while(lines != []):
		for l in lines:
			if l.startswith('#CHROM'):
				g = reduce(lambda x,y: x+','+y, l.strip('\n').split('\t')[9:])
				outputfile.write('\t'+str(g) +'\n')
			if not l.startswith('#'):
				tokens = l.strip('\n').split('\t')
				f = filter(lambda x: 'GP' in x, tokens[7].split(';'))
				if f != []:
					pos = 'chr'+f[0].split('=')[1].split(':')[0]+'pos'+f[0].split('=')[1].split(':')[1]
					ref[pos] = tokens[3]
					alt[pos] = tokens[4]
					m=pos +'\t'
					for t in tokens[9:]:
						m = m + str(int(t[0]) + int(t[2])) + ','
					outputfile.write(m.strip(',')+'\n')
		lines = file.readlines(1000000)
	glob.dump(ref, outputname+'Ref')
	glob.dump(alt, outputname+'Alt')


