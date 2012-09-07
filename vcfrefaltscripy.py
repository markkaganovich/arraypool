import simplejson

outputname ='testceu'
vcffile = '../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf'

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