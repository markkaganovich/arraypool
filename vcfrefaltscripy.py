import simplejson



#find reversed ref/alt genotypes
'''
f = open('./omniexpresssnps')
l= f.readline()
f.close()
omniexpresssnps = map(lambda x: eval(x), l[1:].split(',')[0:-1])

mat = []
reversemat = []
inhapmap = []
inref1kg = []
keys = set(hapmapRefAlt.keys())
ref1kgkeys = set(ref1kg.keys())
for k in omniexpresssnps:
	if k in keys:
		inhapmap.append(k)
	if k in ref1kgkeys:
		inref1kg.append(k)
	try:
		if ref1kg[k] == hapmapRefAlt[k]['ref'] and alt1kg[k] == hapmapRefAlt[k]['alt']:
			mat.append(k)
		elif ref1kg[k] == hapmapRefAlt[k]['alt'] and alt1kg[k] == hapmapRefAlt[k]['ref']:
			reversemat.append(k)
	except KeyError:
		pass

'''
#redo genotypes for those snps that are reverse matches

#one chr at a time
def revgeno(chrom, reversematch):
	SnpPos = map(lambda x: x.strip('\n'), open('hapmapchr'+str(chrom)+'SnpPos').readlines())
	SnpPosset = set(SnpPos)
	geno = map(lambda x: x.strip('\n'), open('hapmapchr'+str(chrom)+'genotype').readlines())
	posi = []
	for s in reversematchset:
		if s in SnpPosset:
			posi.append(SnpPos.index(s))
	posiset = set(posi)
	newgeno = []
	for i in range(0, len(geno)):
		if i in posi:
			g = map(lambda x: int(x), geno[i].split(','))
			newg = map(lambda x: abs(2-x), g)
			ng = ''
			for x in newg:
				ng = ng+str(x)+','
			newgeno.append(ng.strip(','))
		else:
			newgeno.append(geno[i])
	output = open('hapmapchr'+str(chrom)+'genotype','w')
	for n in newgeno:
		output.write(n)
		output.write('\n')
	output.close()	
	
#for c in range(1,23):
#	revgeno(c)

def calcmatches(vcfinputgroup, outputname):
	vcffile = '../1000GenomesData/' + vcfinputgroup + '.low_coverage.2010_09.genotypes.vcf'
	file = open(vcffile)
	refhash = {}
	althash = {}
	lines = file.readlines(1000000)
	outputref = open(outputname + 'Ref','w')
	outputalt = open(outputname + 'Alt','w')
	while(lines != []):
		for l in lines:
			if not l.startswith('#'):
				tokens = l.strip('\n').split('\t')
				f = filter(lambda x: 'GP' in x, tokens[7].split(';'))
				if f != []:
					pos = 'chr'+f[0].split('=')[1].split(':')[0]+'pos'+f[0].split('=')[1].split(':')[1]
					ref = tokens[3]
					alt = tokens[4]
					refhash[pos] = ref
					althash[pos] = alt
		lines = file.readlines(1000000)
	simplejson.dump(refhash, outputref)
	simplejson.dump(althash, outputalt)