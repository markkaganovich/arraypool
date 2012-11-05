from pygr import worldbase
import glob
import os

def filterSNPs(name):
	reffile = name +'Ref'
	altfile = name + 'Alt'
	[ref, alt] = map(lambda x: glob.json(x, ''), [reffile, altfile])
	print "Loaded Ref {0}, Alt {1}".format(len(ref),len(alt))
	keys = ref.keys()
	complsnps = []
	for snppos in keys:
		if glob.compl[ref[snppos].upper()] == alt[snppos].upper() or ref[snppos].upper() == alt[snppos].upper():
			complsnps.append(snppos)
			del ref[snppos]
			del alt[snppos]

	print len(ref)
	print len(alt)
	glob.dump(ref, reffile+'T')
	glob.dump(alt, altfile+'T')
	return complsnps

	
def checkRef(name):
	reffile = name +'RefT'
	altfile = name + 'AltT'
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = True)
	print "Loaded hg19"
	ref = glob.json(reffile)
	print "Loaded Ref"
	alt = glob.json(altfile)
	print "Loaded Alt"
	flip = []
	errors = []
	keys = ref.keys()
	for snppos in keys:
		print snppos + '\t' + name
		t = snppos.split('pos')
		hg19snp = str(hg19[t[0]][int(t[1])-1]).upper()
		refsnp = ref[snppos].upper()
		altsnp = alt[snppos].upper()
		if hg19snp == refsnp:
			continue
		elif hg19snp == altsnp:
			flip.append(snppos)
		else:
			print "Error: Neither Ref nor Alt of SNP corresponds to hg19 sequence"
			errors.append(snppos)
	glob.dump(flip, name+'flips')
	glob.dump(errors, name+'errors')
	return [flip, errors]
	
def corrRef(flip, name):
	reffile = name +'RefT'
	altfile = name + 'AltT'
	for snp in flip:
		t = ref[snp]
		ref[snp] = alt[snp]
		alt[snp] = t
	glob.dump(ref, reffile+'flipped')
	glob.dump(alt, altfile+'flipped')	
		
def flipGeno(genofile, flip):
	"""flip genotypes based on new ref and alt,
	the new genotype allele frequencies are 
	flipped 0->2, 2->0 based on list of snp positions 
	inputed as flip list
	"""	
	
	lines = open(genofile).readlines()
	newgeno = open(genofile+'flipped', 'w')
	newgeno.write(lines[0])
	for l in lines[1:]:
		t = l.split('\t')
		if t[0] in flip:
			print t[1]
			newg = map(lambda x: 2-int(x), t[1].strip('\n').strip(',').split(','))
			newl = t[0] +'\t' 
			newl = reduce(lambda x,y: x+str(y) + ',', [newl]+newg)
			newgeno.write(newl.strip('\n'))
		else:
			newgeno.write(l)
		
def flipArray(arrayname, flip):
	"""flip array snp frequencies (hash)
	1-freq for those in snp list inputed as flip
	input is constructed in the original getarraysnps() function
	""" 
	
	try:
		arrayfreq = glob.load(arrayname+'freq')
	except:
		"No array snp frequency file"
	for snp in flip:
		arrayfreq[flip] = 1 - arrayfreq[flip]
	glob.dump(arrayfreq, arrayname+'freq')
	
def printtabarray(arrayname):
	"""output will be analyzed by R to find cell line frequencies
	"""
	
	output = open(arrayname+'Rinput', 'w')
	freq = glob.load(arrayname+'freq')
	for snp in freq.keys():
		output.write(snp + '\t')
		output.write(str(freq[snp]) + '\n') 
			
#modified to the 19 version
def parse1KGvcf(vcffile, outputname):
	file = open(vcffile)
	outputfile = open(outputname+'Geno', 'w')
	ref = {}
	alt = {}
	lines = file.readlines(1000000)
	while(lines != []):
		for l in lines:
			if l.startswith('#CHROM'):
				g = reduce(lambda x,y: x+','+y, l.strip('\n').split('\t')[9:])
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
	glob.dump(ref, outputname+'Ref')
	glob.dump(alt, outputname+'Alt')
	
def parsehapmap():
	import parsehapmapgenotypes
	ref = {}
	alt = {}	
	genotype = open('hapmapGeno','w')
	for c in range(1,23):
		"here"
		[r,a] = parsehapmapgenotypes.parsehapmapchrom(c)
		ref.update(r)
		alt.update(a)
		with open('../genotypes/hapmapchr'+str(c)+'genotype') as g:
			lines = g.readlines()
			map(lambda l: genotype.write(l), lines)
	genotype.close()
	glob.dump(ref, '../genotypes/hapmapRef')
	glob.dump(alt, '../genotypes/hapmapAlt')
	
# Tests
def test():
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = True)
	for i in range(1,1000000):
		b = str(hg19['chr1'][i])

def testRef():
	reffile = '../genotypes/CEUlowcovRef'
	altfile = '../genotypes/CEUlowcovAlt'
	ref = simplejson.load(open(reffile))
	print "loaded ref"
	alt = simplejson.load(open(altfile))
	print "all loaded"
	keys = ref.keys()
	for snppos in keys:
		continue
	
	