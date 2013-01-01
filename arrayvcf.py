import vcf
import gl
import os
import argparse

'''
the .vcf genotype here is assumed to have an INFO field with hg19 positions encoded, for later vcf files for 1kg
the primary chr pos should be used since it is already hg19

order of implementation: 
		1) run parse1KGvcf to find out Ref / Alt stuff, and get genotype matrix from pool lines and .vcf file
		2) run getarraysnps using above Ref / Alt distinction with
			a) Uniform array 
			b) experiment array
			get final array snp list from this	
		3) use snp list from getarraysnps to remake genotype file for linear regression next step 

'''

def runeverything(uniformarray, exparray, vcffile, pool, refdb, altdb, genotypedb):
	'''
	runeverything('MKReportbySNP1.txt', 'MKReportbySNP3.txt', '../1000GenomesData/low_coverage.merged.vcf', p1lines, 'testoutputRef', 'testoutputAlt', 'testoutput')
		returns all .Rinput files: uniformarray.Rinput, exparray.Rinput, poolgenotype.Rinput

	'''
	poollines = gl.jsonload(pool)
	parse1KGvcf(vcffile, poollines, genotypedb, refdb, altdb)
	arrays(uniformarray, exparray, refdb, altdb, genotypedb)

def arrays(uniformarray, exparray, refdb, altdb, genotypedb, genooutputname):
	'''
	array('MKReportbySNP1.txt', 'MKReportbySNP3.txt', 'testoutputRef', 'testoutputAlt','testoutput')
		returns processed *.Rinput files for each array (control and experiment)
		the input is a genotypedb file extracted from 1kg population .vcf and the refdb altdb from that same .vcf
	'''

	usnplist = getarraysnps(uniformarray, refdb, altdb)
	esnplist = getarraysnps(exparray, refdb, altdb)
	jointsnplist = sorted(list(set(usnplist) & set(esnplist)), key = lambda x: (int(x.split(':')[0]), int(x.split(':')[1])))
	finalsnplist = reshapegenotype(genotypedb, jointsnplist, outputname=genooutputname)
	printtabarray(finalsnplist,uniformarray)
	printtabarray(finalsnplist, exparray)


def parse1KGvcf(vcffile, poollines, genotypedboutput, refdboutput, altdboutput):
	'''
	parse1KGvcf('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf' , p1lines, 'testoutput', 'testoutputRef', 'testoutputAlt')
	
	'''
	vfile = open(vcffile, 'r')
	vcf_reader = vcf.Reader(vfile)
	
	outputfile = open(genotypedboutput, 'w')
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
		m = chrom+':'+pos+'\t'
		sm = 0
		for s in poollines:
			try:				
				geno = record.genotype(s)['GT']
			except KeyError:
				continue
			if geno:
				if '|' in geno:
					g = geno.split('|')
				if '\\' in geno:
				 	g = geno.split('\\')
				if '/' in geno:
					g = geno.split('/')
				an = int(g[0]) + int(g[1])
				sm = sm + an
				m = m + str(an) + ','
			else:
				m = m+str(0) + ','
		if sm > 0:
			outputfile.write(m.strip(',') + '\n')

	gl.jsondump(ref, refdboutput)
	gl.jsondump(alt, altdboutput)
	
	
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
	h = ['SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top', 'GC Score', 'Allele1 - Forward', 'Allele2 - Forward', 'Allele1 - Plus', 'Allele2 - Plus', 'Chr', 'Position', 'GT Score', 'Cluster Sep', 'SNP', 'X', 'Y', 'X Raw', 'Y Raw', 'B Allele Freq', 'Log R Ratio', 'CNV Value', 'CNV Confidence', 'Top Genomic Sequence', 'Plus/Minus Strand', 'Theta', 'R']
	#h = 'SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tGC Score\tSNP Index\tAllele1 - Forward\tAllele2 - Forward\tAllele1 - AB\tAllele2 - AB\tAllele1 - Plus\tAllele2 - Plus\tChr\tPosition\tSNP\tILMN Strand\tTop Genomic Sequence\tPlus/Minus Strand\tTheta\tR\tX\tY\tX Raw\tY Raw\tB Allele Freq\r\n'.split('\t')
	snpi = h.index("SNP")
	chri = h.index("Chr")
	posi = h.index("Position")
	#Yi = h.index("Y")
	#Xi = h.index("X")
	#Ri = h.index("R")
	thetai = h.index("Theta")
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
					#f = float(t[Yi])/(float(t[Ri])) 
					f = float(t[thetai])
				elif genoref[snppos] == alt and genoalt[snppos] == ref:
					#f = 1 - (float(t[Yi])/(float(t[Ri])))
					 f = 1 - float(t[thetai])
				else:
					continue
				if f != 0 and f != 1:
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
 

def reshapegenotype(genofile, arraysnps, outputname = 'poolgenotype3.Rinput'):
	'''
	reshapegenotype('testoutputGeno', 'MKReportbySNP1.txtsnps', 'testgenotype.Rinput')
	
	'''
	genof = open(genofile, 'r')
	arraysnpset = set(arraysnps)
	out = open(outputname, 'w')
	jointlist = []

	for l in genof:
		snppos = l.split('\t')[0]
		if snppos in arraysnpset:
			out.write(snppos+','+l.split('\t')[1])
			jointlist.append(snppos)
	return jointlist

	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-runeverything', action='store_true')
	parser.add_argument('arrays', nargs='+', help="Process arrays: uniform array, experiment array, refdb, altdb, genotypedb, genooutputname, vcffile, poollines")
	parser.add_argument('--parse1KGvcf', action='store_true')
	args = parser.parse_args()

	#if args.runeverything:
	#	runeverything('25M1.1', '25M1.3', '../1000GenomesData/low_coverage.merged.vcf', 'pool1', '25Marrays1230', '25Marrays1230Ref', '25Marrays1230Alt')
	
	controlarray = args.arrays[0]
	print "control array: {0}".format(controlarray)
	exparray = args.arrays[1]
	print "experiment array: {0}".format(exparray)
	refdb = args.arrays[2]
	print "refdb: {0}".format(refdb)
	altdb = args.arrays[3]
	print "altdb: {0}".format(altdb)
	genotypedb = args.arrays[4]
	print "genotypedb: {0}".format(genotypedb)
	genooutputname = args.arrays[5]
	print "outputname for genotype Rinput: {0}".format(genooutputname)

	if args.parse1KGvcf:
		vcffile = args.arrays[6]
		print "vcf file: {0}".format(vcffile)
		poollines = args.arrays[7]
		print "pool lines: {0}".format(poollines)
		parse1KGvcf(vcffile, poollines, genotypedb, refdb, altdb)
	
	arrays(args.arrays[0], args.arrays[1], args.arrays[2], args.arrays[3], args.arrays[4], args.arrays[5])




