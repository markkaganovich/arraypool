import simplejson
import parsegenotypes

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
			
	output = open('omni25Msnpssorted','w')
	simplejson.dump(sortedsnps, output)

### get genotype
class genotypes:
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
	print(len(lines))
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
	file = open(filestruc.name+'tempgenos','w')
	simplejson.dump(genos, file)
	file.close()		
	return genos
			
def combinegenos(names, chosenSNPs):			
	genos = {}
	files = map(lambda x: genotypes(x), names)
	genos = reduce(lambda x,y: getsnpgenos(x, y, chosenSNPs), [genos]+files)
	file = open('Genos','w')
	simplejson.dump(genos, file)
	file.close()
		
names = ['../genotypes/CEUlowcov','../genotypes/YRIlowcov','../genotypes/CHBJPTlowcov','../genotypes/YRItrio', '../genotypes/CEUtrio']		

vcffiles = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf','../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', 
'../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']

""" 
parse genotype files and flip them around according to hg19

parsegenotypes.parse1KGvcf(vcffiles[i], names[i])
c = parsegenotypes.filterSNPs(names[i])
f = parsegenotypes.checkRef(names[i])
parsegenotypes.flipGeno(names[i]+'Geno', f[0])

""" 

file = open('./omni25Msnpssorted')
snps = simplejson.load(file)
file.close()
print len(snps)
combinegenos(names, snps)

