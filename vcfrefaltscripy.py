import simplejson

#find reversed ref/alt genotypes
f = open('../genotypes/omniexpresssnps')
l= f.readline()
f.close()
omniexpresssnps = map(lambda x: eval(x), l[1:].split(',')[0:-1])

#redo genotypes for those snps that are reverse matches

#one chr at a time
def revgeno(chrom, reversematch):
	SnpPos = map(lambda x: x.strip('\n'), open('../genotypes/hapmapchr'+str(chrom)+'SnpPos').readlines())
	SnpPosset = set(SnpPos)
	geno = map(lambda x: x.strip('\n'), open('../genotypes/hapmapchr'+str(chrom)+'genotype').readlines())
	posi = []
	for s in reversematch:
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
	output = open('../genotypes/hapmapchr'+str(chrom)+'genotype','w')
	for n in newgeno:
		output.write(n)
		output.write('\n')
	output.close()	
	

def calcmatches(vcfinputgroup, outputname):
	vcffile = '../1000GenomesData/' + vcfinputgroup + '.2010_09.genotypes.vcf'
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
	
# test if ceu and yri ref/alt designations are the same:
'''
f = open('../genotypes/testceuRef')
ref1kg = simplejson.load(f)
f.close()
f = open('../genotypes/testceuAlt')
alt1kg = simplejson.load(f)
f.close()
ref2kg = simplejson.load(open('../genotypes/YRIRef'))
alt2kg = simplejson.load(open('../genotypes/YRIAlt'))

'''
def calcreversals(omniexpresssnps, ref1, alt1, ref2, alt2):
	Ref = {}
	Alt = {}
	mat = []
	reversemat = []
	inhapmap = []
	inref1kg = []
	keys1 = set(ref1.keys())
	keys2 = set(ref2.keys())
	for k in omniexpresssnps:
		if k in keys1:
			inhapmap.append(k)
		if k in keys2:
			inref1kg.append(k)
		try:
			if ref1[k] == ref2[k] and alt1[k] == alt2[k]:
				mat.append(k)
				Ref[k] = ref1[k]
				Alt[k] = alt1[k]
			elif ref1[k] == alt2[k] and alt1[k] == ref2[k]:
				reversemat.append(k)
				Ref[k] = ref1[k]
				Alt[k] = alt1[k]
		except KeyError:
			pass
	return [Ref, Alt, mat, reversemat]
	
#r = calcreversals(omniexpresssnps, ref1kg, alt1kg, ref2kg, alt2kg)

#Ref = {}
#Alt = {}	
def mdic(n):
	# type can be 'Ref' or 'Alt'
	# this function opens dictionaries
	print n
	f = open(n+'Alt')
	ref = simplejson.load(f)
	f.close()
	return ref
#names = ['testceu', 'YRI', 'CHBJPT', 'CEUtrio', 'YRItrio']
#l = map(mdic, names)
#for n in l:
#	Ref.update(n)
# save as 1kgRef and 1kgAlt

# add hapmap to Ref and Alt and flip genotypes for each chrom
# aggregate all hapmap chroms

'''
Ref = {}
Alt = {}
for c in range(1,23):
	f = open('../genotypes/hapmapchr'+str(c)+'RefAlt')
	refaltdic = simplejson.load(f)
	f.close()
	keys = refaltdic.keys()
	for k in keys:
		Ref[k] = refaltdic[k]['ref']
		Alt[k]  = refaltdic[k]['alt']
'''
#r = calcreversals(omniexpresssnps, ref1kg, alt1kg, hapmapRef, hapmapAlt)
	
#for c in range(1,23):
#	revgeno(c, reversematches)	
		
	
	

