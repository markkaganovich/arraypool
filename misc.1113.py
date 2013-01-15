'''
misc code for arrayvcf.py that wasn't really used

'''
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

def comparearrays(array1, array2):
	thetai = 10
	snpidi = 0
	a1 = open(array1)
	a2 = open(array2)
	out = open(array1+array2 +'.diff.theta.R', 'w')
	for l1, l2 in zip(a1, a2):
			t1 = l1.split('\t')
			t2 = l2.split('\t')
			if t1[snpidi] == t2[snpidi] and t1[thetai] != 0:
				try:
					out.write(str((float(t1[thetai]) - float(t2[thetai]))/sum([float(t1[thetai]), float(t2[thetai])])/2) +'\n')
				except ZeroDivisionError:
					continue
	out.close()



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




'''
	parser = argparse.ArgumentParser()
	parser.add_argument('arrays', nargs='+', help="Process arrays: arrays, refdb, altdb, output_ext")
	parser.add_argument('-chr')

	parser.add_argument('--parse1KGvcf', action='store_true', help="Make genotypedb matrix: vcffile, poollines, genotpedbname")
	args = parser.parse_args()

	print args
	
	if args.parse1KGvcf:
		vcffile = args.arrays[0]
		print "vcf file: {0}".format(vcffile)
		poollines = args.arrays[1]
		print "pool lines: {0}".format(poollines)
		genotypedb = args.arrays[2]
		parse1KGvcf(vcffile, poollines, genotypedb, genotypedb+'Ref', genotypedb+'Alt')
	else:
		array = args.arrays[0]
		print "control array: {0}".format(array)
		refdb = args.arrays[1]
		print "refdb: {0}".format(refdb)
		altdb = args.arrays[2]
		print "altdb: {0}".format(altdb)
		genotypedb = args.arrays[3]
		print "output_ext: {0}".format(genotypedb)
		if args.chr:


		arrays(args.arrays[0], args.arrays[1], args.arrays[2], args.arrays[3], chr=args.chr)

'''	