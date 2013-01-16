import vcf
import gl
import os
import argparse
import sys

'''
the .vcf genotype here is assumed to have an INFO field with hg19 positions encoded, for later vcf files for 1kg
the primary chr pos should be used since it is already hg19

order of implementation: 
		1) run parse1KGvcf to find out Ref / Alt stuff, and get genotype matrix from pool lines and .vcf file
		2) run getarraysnps using above Ref / Alt distinction with inputed arrays
		3) in R: use snp list from getarraysnps to remake genotype file for linear regression next step 

'''
def arrays(*args):
	'''
	array('MKReportbySNP1.txt', 'MKReportbySNP3.txt', 'testoutputRef', 'testoutputAlt','testoutput')
		returns processed *.Rinput files for each array (control and experiment)
		the input is a genotypedb file extracted from 1kg population .vcf and the refdb altdb from that same .vcf
	'''
	array = args[0]
	refdb = args[1]
	altdb = args[2]
	output_ext = args[3]
	kwargs = args[4]
	getarraysnps(array, refdb, altdb, array+output_ext, kwargs)


def parse1KGvcf(vcffile, plines, genotypedboutput, refdboutput, altdboutput):
	'''
	parse1KGvcf('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf' , p1lines, 'testoutput', 'testoutputRef', 'testoutputAlt')
	'''

	vfile = open(vcffile, 'r')
	vcf_reader = vcf.Reader(vfile)
	
	poollines = gl.jsonload(plines)
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
			record.genotype(s)
			try:				
				geno = record.genotype(s)['GT']
			except KeyError:
				print "no genotype"
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
	
	
def getarraysnps(report, fgenoref, fgenoalt, outname, kwargs):
	'''
	test this with MKReport1bysnps.txt and output of parse1KGvcf function
	getarraysnps('MKReport1bysnps.txt', 'testoutputRef', 'testoutputAlt')
	In this version we totally ignore compliments

	kwargs:
	snp
	chr
	pos
	theta
	header
	'''

	print report
	print kwargs
	file = open(report)
	lines = file.readlines()
	file.close()

	genoref = gl.jsonload(fgenoref)  
	genoalt = gl.jsonload(fgenoalt)

	output = open(outname, 'w')

	try:
		header = kwargs['header']
		h = header.split('\\t')
		print h
		snpi = h.index("SNP")
		chri = h.index("Chr")
		posi = h.index("Position")
		thetai = h.index("Theta")
	except KeyError:
		try:
			snpi = kwargs['snp']
			chri = kwargs['chr']
			posi = kwargs['pos']
			thetai = kwargs['theta']
		except KeyError:
			print "No header or column numbers provided provided"

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
				if genoref[snppos] == ref and genoalt[snppos] == alt:
					f = float(t[thetai])
				elif genoref[snppos] == alt and genoalt[snppos] == ref:
					 f = 1 - float(t[thetai])
				else:
					continue
				if f != 0 and f != 1:
					freq[snppos] = f
					snplist.append(snppos)
					output.write(snppos + '\t' + str(f) + '\n')
		except:
			continue
		
	gl.jsondump(snplist, report+'snps')
	gl.jsondump(freq, report+'freq')

	
if __name__ == "__main__":
	
	args = sys.argv[1:]
	print args 
	if '-parse1KGvcf' in args:
		vcffile = args[1]
		poollines = args[2]
		genotypedb = args[3]
		parse1KGvcf(vcffile, poollines, genotypedb, genotypedb+'Ref', genotypedb+'Alt')
	if 'arrays' in args:
		#array = args[1]
		#refdb = args[2]
		#altdb = args[3]
		#output_ext = args[4]
		kwargs = args[5:]
		k = dict(map(lambda x: x.split('='), kwargs))

		#arrays(array, refdb, altdb, output_ext, k)
		arrays(args[1], args[2], args[3], args[4], k)





