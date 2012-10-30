import simplejson
import parsegenotypes
import globals

def getarraysnps():
	report = 'Array25M1'
	file = open(report)
	lines = file.readlines()
	file.close()
	header = "SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,Allele1 - Plus,Allele2 - Plus,Chr,Position,SNP,Theta,R,X,Y,X Raw,Y Raw,B Allele Freq"
	t = header.split(',')
	snpi = t.index("SNP")
	chri = t.index("Chr")
	posi = t.index("Position")
	Yi = t.index("Y")
	Ri = t.index("R")
	snplist = []
	for l in lines:
		t = l.split('\t')
		snplist.append('chr'+t[chri]+'pos'+t[posi])
	#sortedsnpschr = sorted(snplist, key=lambda snp: snp.split('pos')[0].split('chr')[1])
	sortedsnps = []
	for c in range(1,23):
		temp = []
		for s in snplist:
			if str(c) == s.split('pos')[0].split('chr')[1]:
				temp.append(s)
		sortedsnps.extend(sorted(temp, key=lambda snp: int(snp.split('pos')[1])))
	globals.dump(sortedsnps, 'omni25Msnpssorted')

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
		
names = ['../genotypes/CEUlowcov','../genotypes/YRIlowcov','../genotypes/CHBJPTlowcov','../genotypes/YRItrio', '../genotypes/CEUtrio']		

vcffiles = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf','../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', 
'../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']

""" 
#parse genotype files and flip them around according to hg19
for i in range(0, len(names)):
	parsegenotypes.parse1KGvcf(vcffiles[i], names[i])
	b = parsegenotypes.filterSNPs(names[i])
	print "{0} SNPs filtered out".format(len(b))
	c = parsegenotypes.checkRef(names[i])
	print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
	parsegenotypes.flipGeno(names[i]+'Geno', c[0])



snps = globals.json('omni25Msnpssorted')
print len(snps)
combinegenos(names, snps)
"""
#b = parsegenotypes.filterSNPs('../genotypes/hapmap')
#print "{0} SNPs filtered out".format(len(b))
c = parsegenotypes.checkRef('../genotypes/hapmap')
print "{0} errors and {1} flipped".format(len(c[1]),len(c[0]))
parsegenotypes.flipGeno('../genotypes/hapmapgenotype', c[0])

