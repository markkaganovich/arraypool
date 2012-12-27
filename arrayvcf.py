import vcf
import gl
import os
import argparse

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
'''
order: 1) run parse1KGvcf to find out Ref / Alt stuff, and get genotype matrix from pool lines and .vcf file
		2) run getarraysnps using above Ref / Alt distinction with
			a) Uniform array 
			b) experiment array
			get final array snp list from this	
		3) use snp list from getarraysnps to remake genotype file for linear regression next step 

'''

def runeverything(uniformarray, exparray, vcffile, pool, output):
	'''
	runeverything('MKReportbySNP1.txt', 'MKReportbySNP3.txt', '../1000GenomesData/low_coverage.merged.vcf', p1lines, 'testoutput')
		returns all .Rinput files: uniformarray.Rinput, exparray.Rinput, poolgenotype.Rinput
vc
	'''
	poollines = gl.jsonload(pool)
	parse1KGvcf(vcffile , output, poollines)
	usnplist = getarraysnps(uniformarray, 'testoutputRef', 'testoutputAlt')
	esnplist = getarraysnps(exparray, 'testoutputRef', 'testoutputAlt')
	
	jointsnplist = list(set(usnplist) & set(esnplist))

	printtabarray(jointsnplist,uniformarray)
	printtabarray(jointsnplist, exparray)

	reshapegenotype(output+'Geno', jointsnplist)



def parse1KGvcf(vcffile, outputname, poollines):
	'''
	parse1KGvcf('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf' , 'testoutput', p1lines)
	
	'''
	vfile = open(vcffile, 'r')
	vcf_reader = vcf.Reader(vfile)
	
	outputfile = open(outputname+'Geno', 'w')
	ref = {}
	alt = {}
	
	for record in vcf_reader:
		try:
			chrom = record.INFO['GP'].split(':')[0]
			pos = record.INFO['GP'].split(':')[1]	
		except KeyError:
			continue
		ref[chrom+':'+pos] = str(record.REF) 
		alt[chrom+':'+pos] = str(record.ALT[0])
		poolsamples = filter(lambda x: x.sample in poollines, record.samples)
		m = chrom+':'+pos+'\t'
		for s in poolsamples:
			if s['GT']:
				if '|' in s['GT']:
					g = s['GT'].split('|')
				if '\\' in s['GT']:
				 	g = s['GT'].split('\\')
				if '/' in s['GT']:
					g = s['GT'].split('/')
				m = m + str(int(g[0]) + int(g[1])) + ','
			else:
				m = m+str(0) + ','
		outputfile.write(m.strip(',') + '\n')

	gl.jsondump(ref, outputname+'Ref')
	gl.jsondump(alt, outputname+'Alt')
	
	#checkRef
	
def getarraysnps(report, fgenoref, fgenoalt):
	'''
	test this with MKReport1bysnps.txt and output of parse1KGvcf function

	getarraysnps('MKReport1bysnps.txt', 'testoutputRef', 'testoutputAlt')

	In this version we totally ignore compliments
	'''

	print report
	file = open(report)
	lines = file.readlines()
	file.close()

	genoref = gl.jsonload(fgenoref)  # later can make this in-memory
	genoalt = gl.jsonload(fgenoalt)
	#header = "SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,Allele1 - Plus,Allele2 - Plus,Chr,Position,SNP,Theta,R,X,Y,X Raw,Y Raw,B Allele Freq"
	#header = "SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,Allele1 - Forward,Allele2 - Forward,Allele1 - Plus,Allele2 - Plus,Chr,Position,GT Score,Cluster Sep,SNP,X,Y,X Raw,Y Raw,B Allele,Freq,Log R Ratio,CNV Value,CNV Confidence"
	#h = header.split(',')
	#h = ['SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top', 'GC Score', 'Allele1 - Forward', 'Allele2 - Forward', 'Allele1 - Plus', 'Allele2 - Plus', 'Chr', 'Position', 'GT Score', 'Cluster Sep', 'SNP', 'X', 'Y', 'X Raw', 'Y Raw', 'B Allele Freq', 'Log R Ratio', 'CNV Value', 'CNV Confidence', 'Top Genomic Sequence', 'Plus/Minus Strand', 'Theta', 'R\r\n']
	h = 'SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tGC Score\tSNP Index\tAllele1 - Forward\tAllele2 - Forward\tAllele1 - AB\tAllele2 - AB\tAllele1 - Plus\tAllele2 - Plus\tChr\tPosition\tSNP\tILMN Strand\tTop Genomic Sequence\tPlus/Minus Strand\tTheta\tR\tX\tY\tX Raw\tY Raw\tB Allele Freq\r\n'.split('\t')
	snpi = h.index("SNP")
	chri = h.index("Chr")
	posi = h.index("Position")
	Yi = h.index("Y")
	Xi = h.index("X")
	snplist = []
	freq = {}
	for l in lines:
		t = l.split('\t')
		try:
			if t[chri] not in map(lambda x: str(x), range(1,23)):
				continue
			else:
				snppos = t[chri]+':'+t[posi]
				ref = t[snpi].split('/')[0][1] 
				alt = t[snpi].split('/')[1][0]
				#f = 0
				if genoref[snppos] == ref and genoalt[snppos] == alt:
					f = float(t[Yi])/(float(t[Yi])+float(t[Xi])) 
				elif genoref[snppos] == alt and genoalt[snppos] == ref:
					f = 1 - (float(t[Yi])/(float(t[Yi])+float(t[Xi])))
				else:
					continue
				if f != 0:
					freq[snppos] = f
					snplist.append(snppos)
		except:
			continue
		
	gl.jsondump(snplist, report+'snps')
	gl.jsondump(freq, report+'freq')

	return snplist

def printtabarray(jointsnplist, arrayname):
		"""output will be analyzed by R to find cell line frequencies
		"""	
		output = open(arrayname+'.Rinput', 'w')
		freq = gl.jsonload(arrayname+'freq')
		for snp in jointsnplist:
			output.write(snp + '\t')
			output.write(str(freq[snp]) + '\n') 
 

def reshapegenotype(genofile, arraysnps, outputname = 'poolgenotype.Rinput'):
	'''
	reshapegenotype('testoutputGeno', 'MKReportbySNP1.txtsnps', 'testgenotype.Rinput')
	
	'''
	genof = open(genofile, 'r')
	arraysnpset = set(arraysnps)
	out = open(outputname, 'w')

	for l in genof:
		snppos = l.split('\t')[0]
		if snppos in arraysnpset:
			out.write(l)



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

	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-runeverything', action='store_true')
	args = parser.parse_args()

	if args.runeverything:
		runeverything('MKReportbySNP1.txt', 'MKReportbySNP3.txt', '../1000GenomesData/low_coverage.merged.vcf', 'pool1', 'testoutput')



