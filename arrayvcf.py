import vcf
import gl
import os

#poollines = ['NA19211', 'NA18943', 'NA19209', 'NA18526' ,'NA19161' , 'NA11920', 'NA11995' , 'NA18564' , 'NA18499' , 'NA12003']
poollines = ["NA18516", "NA18517", "NA18579", "NA18592", "NA18561", "NA07357", "NA06994", "NA18526", "NA12004", "NA19141", "NA19143", "NA19147", "NA19152", "NA19153", "NA19159", "NA19171", "NA19172", "NA19190", "NA19207", "NA19209", "NA19210", "NA19225", "NA19238", "NA19239", "NA18856", "NA18858", "NA18562", "NA18563", "NA18853", "NA18861"]
names1KG = ['CEUlowcov','YRIlowcov','CHBJPTlowcov','YRItrio', 'CEUtrio']	
vcffiles = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf','../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', 
'../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']


def in1kg(poollines, names):
	lines = []
	for n in names:
		lines.extend(getlines(n))
	result = []
	for l in poollines:
		if l in lines:
			result.append(True)
		else:
			result.append(False)
	return result

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
	
		#checkRef
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


def splitreport(f, dir):
	r = open(f)
	header = None
	for l in r:
		if "[Data]" in l:
			header = r.next().split('\t')
			print header
			continue
		if header:
			i = header.index('Sample ID')
			outputfile = '25M' + l.split('\t')[i]
			if outputfile in os.listdir(dir):
				out.write(l)
			else:
				out = open(outputfile, 'w')
				out.write(l)
				print out

	




