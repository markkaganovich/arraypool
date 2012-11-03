import simplejson
import parsegenotypes
import globals

def getarraysnps():
	report = 'Array25M1'
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
		snppos = 'chr'+t[chri]+'pos'+t[posi]
		snplist.append(snppos)
		ref[snppos] = t[snpi].split('/')[0][1] 
		alt[snppos] = t[snpi].split('/')[1][0]
		freq[snppos] = float(t[Yi])/float(t[Ri])
		
	#sortedsnpschr = sorted(snplist, key=lambda snp: snp.split('pos')[0].split('chr')[1])
	"""
	don't think they need to be sorted
	
	sortedsnps = []
	for c in range(1,23):
		temp = []
		for s in snplist:
			if str(c) == s.split('pos')[0].split('chr')[1]:
				temp.append(s)
		sortedsnps.extend(sorted(temp, key=lambda snp: int(snp.split('pos')[1])))
	globals.dump(sortedsnps, 'omni25Msnpssorted')
	"""
	globals.dump(snplist, report+'snps')
	globals.dump(ref, report+'ref')
	globals.dump(alt, report+'alt')
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
		
homedir = '/srv/gs1/projects/snyder/mark'
names = ['/genotypes/CEUlowcov','/genotypes/YRIlowcov','/genotypes/CHBJPTlowcov','/genotypes/YRItrio', '/genotypes/CEUtrio']	
homedirnames = map(lambda x: homedir+x, names)	

vcffiles = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf','../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', 
'../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']


#parse genotype files and flip them around according to hg19
for i in range(0, len(names)):
	#parsegenotypes.parse1KGvcf(vcffiles[i], names[i])
	#b = parsegenotypes.filterSNPs(names[i])
	#print "{0} SNPs filtered out".format(len(b))
	c = parsegenotypes.checkRef(names[i])
	print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
	parsegenotypes.corrRef(c[0], names[i])
	parsegenotypes.flipGeno(names[i]+'Geno', c[0])
"""


snps = globals.json('omni25Msnpssorted')
print len(snps)
combinegenos(names, snps)
"""
#b = parsegenotypes.filterSNPs('../genotypes/hapmap')
#print "{0} SNPs filtered out".format(len(b))
#c = parsegenotypes.checkRef('../genotypes/hapmap')
#print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
#c = globals.json('../genotypes/hapmapflips')
#parsegenotypes.flipGeno('../genotypes/hapmapgenotype', c[0])

