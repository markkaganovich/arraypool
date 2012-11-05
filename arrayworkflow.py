import parsegenotypes
import sys, getopt
import globals

homedir = '/srv/gs1/projects/snyder/mark/genotypes/'
names = ['CEUlowcov','YRIlowcov','CHBJPTlowcov','YRItrio', 'CEUtrio']	
homedirnames = map(lambda x: homedir+x, names)	

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
		
	globals.dump(snplist, report+'snps')
	globals.dump(ref, report+'RefT')
	globals.dump(alt, report+'AltT')
	globals.dump(freq, report+'freq')
	
### get genotype
class Genotypes:
	def __init__(self, name):
		self.name = name
		self.genofile = open('./' + name + 'Geno')
		self.__genolen__()
	def __genolen__(self):
		l = self.genofile.readline()
		self.ln = len(l.split(','))
		self.genofile.seek(0)

def getsnpgenos(genos, filestruc, chosenSNPs):
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
	for s in notingeno:
		try:
			genos[s] = genos[s].strip(',') + ','+('0,' * filestruc.ln)
		except KeyError:
			genos[s] = '0,' * filestruc.ln
	globals.dump(genos, 'tempgenos')				
	return genos
			
def combinegenos(names, chosenSNPs):			
	genos = {}
	files = map(lambda x: Genotypes(x), names)
	genos = reduce(lambda x,y: getsnpgenos(x, y, chosenSNPs), [genos]+files)
	globals.dump(genos, 'Genos')
		
def processhapmap():
	parsegenotypes.parsehapmap()
	b = parsegenotypes.filterSNPs('../genotypes/hapmap')
	print "{0} SNPs filtered out".format(len(b))
	c = parsegenotypes.checkRef('../genotypes/hapmap')
	print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
	#flips = globals.json('../genotypes/hapmapflips')
	parsegenotypes.flipGeno('hapmapGeno', c[0])

def processgenotypes():
	"""parse genotype files and flip them around according to hg19
	
	"""
	for i, n in enumerate(homedirnames):
		parsegenotypes.parse1KGvcf(vcffiles[i], names[i])
		b = parsegenotypes.filterSNPs(names[i])
		print "{0} SNPs filtered out".format(len(b))
		c = parsegenotypes.checkRef(n)
		print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
		#parsegenotypes.corrRef(c[0], n)
		parsegenotypes.flipGeno(n+'Geno', c[0])
			
def main(argv):
	opts, args = getopt.getopt(argv,"a:gh",["report="])
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
			genofiles = names.append('hapmap')
			snps = globals.json('Array251Msnps')
			combinegenos(genofiles, snps)
			
if __name__ == "__main__":
	args = sys.argv[1:]
	if len(args) > 0:
		main(args)
	

