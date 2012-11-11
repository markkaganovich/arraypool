import parsegenotypes
import sys, getopt
import argparse
import glob

homedir = '/srv/gs1/projects/snyder/mark/genotypes/'
names1KG = ['CEUlowcov','YRIlowcov','CHBJPTlowcov','YRItrio', 'CEUtrio']	
homedirnames = map(lambda x: homedir+x, names1KG)	

vcffiles = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf','../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', 
'../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']

def getarraysnps(report):
	print report
	file = open(report)
	lines = file.readlines()
	file.close()
	header = "SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,Allele1 - Plus,Allele2 - Plus,Chr,Position,SNP,Theta,R,X,Y,X Raw,Y Raw,B Allele Freq"
	h = header.split(',')
	snpi = h.index("SNP")
	chri = h.index("Chr")
	posi = h.index("Position")
	Yi = h.index("Y")
	Ri = h.index("R")
	snplist = []
	ref = {}
	alt = {}
	freq = {}
	for l in lines:
		t = l.split('\t')
		if t[chri] not in map(lambda x: str(x), range(1,23)):
			continue
		else:
			snppos = 'chr'+t[chri]+'pos'+t[posi]
			snplist.append(snppos)
			ref[snppos] = t[snpi].split('/')[0][1] 
			alt[snppos] = t[snpi].split('/')[1][0]
			freq[snppos] = float(t[Yi])/float(t[Ri])
		
	glob.dump(snplist, report+'snps')
	glob.dump(ref, report+'RefT')
	glob.dump(alt, report+'AltT')
	glob.dump(freq, report+'freq')
	
### get genotype
class Genotypes:
	def __init__(self, name, ext = "Geno"):
		self.name = name
		self.genofile = open('./' + name + ext)
		self.__genolen__()
	def __genolen__(self):
		l = self.genofile.readline()
		self.ln = len(l.split(','))
		self.genofile.seek(0)

def getsnpgenos(genos, filestruc, chosenSNPs, incarray = 0):
	lines = filestruc.genofile.readlines()
	snppos = map(lambda x: x.split('\t')[0], lines[1:])
	inboth = set(snppos) & set(chosenSNPs)
	notingeno = set(filter(lambda x: x not in inboth, chosenSNPs))
	try:
		genos['lines'] = genos['lines'] +lines[0].split('\t')[1].split(',')
	except KeyError:
		genos['lines'] = lines[0].split('\t')[1].split(',')
	print("Number of Lines in Genotype:" + str(len(lines)))
	for l in lines[1:]:
		t = l.split('\t')
		snp = t[0]
		if snp in inboth:
			try:
				genos[snp] = genos[snp].strip(',') + ','+t[1].strip('\n')
			except KeyError:
				genos[snp] = t[1].strip('\n')
	if incarray == 0:
		for s in notingeno:
			try:
				genos[s] = genos[s].strip(',') + ','+('0,' * filestruc.ln)
			except KeyError:
				genos[s] = '0,' * filestruc.ln
	glob.dump(genos, 'tempgenos')				
	return genos
			
def combinegenos(names, chosenSNPs, out = 'combGenosfile', incarray = 0):			
	genos = {}
	
	if type(names) is str:
		f = Genotypes(names)
		genos = getsnpgenos(genos, f, chosenSNPs, incarray)

	if type(names) is list:
		files = map(lambda x: Genotypes(x), names)
		genos = reduce(lambda x,y: getsnpgenos(x, y, chosenSNPs, incarray), [genos]+files)
	
	glob.dump(genos, out+'.json')
		
	genos['lines'] = map(lambda x: x.strip('\n'), genos['lines'])
	output = open(out, 'w')
	linenames = reduce(lambda x,y: x +',' + y, genos['lines'])
	output.write('\t' + linenames +'\n')
	for g in genos.keys():
		if g != 'lines':
			output.write(g + '\t' + genos[g].strip(',') + '\n')				
			
def interset(genotypes):
	"""the combine genos stuff is predicated on the idea that if a SNP isn't mentioned
	it's genotype is 0
	in reality, if the snps are genotypes as part of different projects (like 1KG vs. HapMap)
	then if a SNP is not present it may be because it was never tested

	This functions finds the intersecting set of snps between two genotype files, which can then be
	used as part of the "chosen SNPs" to select the ultimate final combined Genos file
	"""
	snplist = map(lambda x: getsnps(x), genotypes)
	print len(snplist)
	ineverything = reduce(lambda x,y: set(x) & set(y), snplist)
	return ineverything

def getsnps(genotypefile):
	f = open(genotypefile)
	lines = f.readlines()
	print len(lines)
	f.close()
	snppos = map(lambda x: x.split('\t')[0], lines[1:])
	return snppos			
			
def intercomb(genotypes, out = 'intercomb'):
	"""the purpose of this function is to do a simple genotype file combine, based on select snps, the
	assumption is that the files are from two different categories of genotyping projects, so the snps
	used are an output of the interset() function
	"""
	snps = interset(genotypes)
	combinedgenos = {}
	combinedgenos['lines'] = []
	for g in genotypes:
		lines = open(g).readlines()
		combinedgenos['lines'].append(lines[0].split('\t')[1])
		for l in lines[1:]:
			t = l.split('\t')
			snp = t[0]	
			if snp in snps:
				try:
					combinedgenos[snp] = combinedgenos[snp].strip(',') + ','+t[1].strip('\n')
				except KeyError:
					combinedgenos[snp] = t[1].strip('\n')
		
	output = open(out, 'w')
	combinedgenos['lines'] = map(lambda x: x.strip('\n'), combinedgenos['lines'])
	linenames = reduce(lambda x,y: x +',' + y, combinedgenos['lines'])
	output.write('\t' +linenames + '\n')
	for g in combinedgenos.keys():
		if g != 'lines':
			output.write(g + '\t' + combinedgenos[g].strip(',') + '\n')
	
		
def processhapmap():
	print "parsing hapmap genotypes chrom files"
	parsegenotypes.parsehapmap()
	print "now filtering SNPs"
	b = parsegenotypes.filterSNPs('../genotypes/hapmap')
	print "{0} SNPs filtered out".format(len(b))
	[flips, errors] = parsegenotypes.checkRef('../genotypes/hapmap')
	print "{0} flips and {1} errors".format(len(flips),len(errors))
	#flips = glob.json('../genotypes/hapmapflips')
	#errors = glob.json('../genotypes/hapmaperrors')
	parsegenotypes.flipGeno('hapmapGeno', flips, errors)

def processgenotypes():
	"""parse genotype files and flip them around according to hg19
	
	"""
	for i, n in enumerate(homedirnames):
		parsegenotypes.parse1KGvcf(vcffiles[i], names1KG[i])
		b = parsegenotypes.filterSNPs(names1KG[i])
		print "{0} SNPs filtered out".format(len(b))
		c = parsegenotypes.checkRef(n)
		print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
		#parsegenotypes.corrRef(c[0], n)
		parsegenotypes.flipGeno(n+'Geno', c[0], c[1])
		
def getpoollines(genofile, pool, out = "poolgenotype"):
	output = open(out, 'w')
	
	pool = glob.json(pool)
	g = open(genofile)
	lines = g.readlines()
	g.close()
	
	linenames = lines[0]
	ln = linenames.split('\t')[1].split(',')
	ln = map(lambda x: x.strip('\n'), ln)
	print pool
	poolinds = map(lambda x: ln.index(x), pool)
	
	nl = [ln[i] for i in poolinds]
	newlinenames = ','
	newlinenames = reduce(lambda x,y: x + ',' + y, nl)
	output.write(newlinenames + '\n')
	
	for l in lines[1:]:
		t = l.split('\t')
		
		gs = t[1].split(',')
		gs = map(lambda x: x.strip('\n'), gs)
		if len(filter(lambda x: int(x)!=0, gs)) >0:
			newl = t[0] + ','
			newg = []
			for i in range(0, len(gs)):
				if i in poolinds:
					newg.append(gs[i])
			newg = reduce(lambda x,y: x + ',' + y, newg)
			newl = newl + newg + '\n'
			output.write(newl)		
			
def mergearraypool(poolgenotypefile, arrayfreqfile):
	pfile = open(poolgenotypefile)
	plines = pfile.readlines()
	pfile.close()
	afile = open(arrayfreqfile)
	alines = afile.readlines()
	afile.close()
	asnps = map(lambda x: x.split('\t')[0], alines[1:])
	asnpsset = set(asnps)
	
	print asnpsset
	output = open("mergedarraypool.Rinput",'w')
	for g in plines[1:]:
		print g.split(',')[0]
		if g.split(',')[0] in asnpsset:
			output.write(g)	
			
			
def main(argv):
	opts, args = getopt.getopt(argv,"a:ghc:i",["report=", "genonames="])
	for opt, arg in opts:
		if opt == '-a':
			report = arg
			print "processing array, getting snps {0}".format(report)
			getarraysnps(arg)
			[f, e] = parsegenotypes.checkRef(report)
			parsegenotypes.flipArray(report, f)
			parsegenotypes.filterzeros(report)
			parsegenotypes.printtabarray(report)
			
		elif opt == '-g':
			processgenotypes()
		elif opt == '-h':
			processhapmap()
		elif opt == '-c':
			snps = glob.json('Array251Msnps')
			if arg == '':
				genofiles = names1KG
			else:
				genofiles = arg
			combinegenos(genofiles, snps)
			
		elif opt == "-i":
			genofiles = names.append('hapmap')
			
			
if __name__ == "__main__":
	#args = sys.argv[1:]
	#if len(args) > 0:
	#	main(args)
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', help='Extract snps from this Array')
	parser.add_argument('-g', action='store_true')
	parser.add_argument('-hapmap', action='store_true')
	parser.add_argument('-init1KG', action='store_true')
	parser.add_argument('-inithapmap', action='store_true')
	parser.add_argument('-int', action='store_true')
	parser.add_argument('-pool', action='store_true')
	parser.add_argument('-i')
	args = parser.parse_args()
	print args
	
	if args.a:
		print "processing array, getting snps {0}".format(args.a)
		report = args.a
		getarraysnps(report)
		[f, e] = parsegenotypes.checkRef(report)
		parsegenotypes.flipArray(report, f, e)
		parsegenotypes.filterzeros(report)
		parsegenotypes.printtabarray(report)
	if args.g:
		processgenotypes()
	if args.hapmap:
		processhapmap()
	if args.init1KG:
		snps = glob.json('Array25M1snps')
		combinegenos(names1KG, snps, 'Genos1kgArray25M')
	if args.inithapmap:
		snps = glob.json('Array25M1snps')
		combinegenos('hapmap', snps, 'hapmapGenosArray25M', 1)
		
	if args.pool:
		getpoollines('intercomb','pool2', 'intercombpool2')
	
	if args.int:
		intercomb(['Genos1kgArray25M', 'hapmapGenosArray25M'])		
		
	
	
	
	

