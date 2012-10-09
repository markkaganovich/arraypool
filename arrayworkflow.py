import simplejson

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

def getsnpgenos(filestruc, chosenSNPs, genos = {}):
	lines = filestruc.genofile.readlines()
	genos = {}
	snppos = map(lambda x: x.split('\t')[0], lines[1:])
	inboth = set(snppos) & set(chosenSNPs)
	notingeno = set([x if x not in inboth for x in chosenSNPs])
	for l in lines[1:]:
		t = l.split('\t')
		snp = t[0]
		if snp in inboth:
			try:
				genos[snp] = genos[snp] + t[1].strip('\n')
			except KeyError:
				genos[snp] = t[1].strip('\n')
		elif snp in notingeno:
			try:
				genos[snp] = genos[snp] + '0,' * filestruc.ln
			except KeyError:
				genos[snp] = '0,' * filestruc.ln
			
def combinegenos()			
	genos = {}
	files = map(lambda x: genotypes(x), names)
	map(lambda x: getsnpgenos(x, chosenSNPs, genos), files)
	file = open('tempdic','w')
	simplejson.dump(genos, file)
	file.close()
			
		
		
		
	
	
#modifyied to the 19 version
def flatfilevcf(vcffile, outputname):
	file = open(vcffile)
	outputfile = open(outputname+'Geno', 'w')
	ref = {}
	alt = {}
	lines = file.readlines(1000000)
	outputref = open(outputname + 'Ref','w')
	outputalt = open(outputname + 'Alt', 'w')
	while(lines != []):
		for l in lines:
			if l.startswith('#CHROM'):
				genomenames = l.strip('\n').split('\t')
				g = [',' if x == '\t' else x for x in genomenames[9:]]
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
	simplejson.dump(ref, outputref)
	simplejson.dump(alt, outputalt)
