from pygr import worldbase
import globals

def filterSNPs(name):
	reffile = name +'Ref'
	altfile = name + 'Alt'
	[ref, alt] = map(lambda x: globals.json(x, ''), [reffile, altfile])
	print "Loaded Ref, Alt"
	print len(ref)
	print len(alt)
	keys = ref.keys()
	complsnps = []
	for snppos in keys:
		if globals.compl[ref[snppos].upper()] == alt[snppos].upper():
			complsnps.append(snppos)
			del ref[snppos]
			del alt[snppos]
	print len(ref)
	print len(alt)
	map(lambda x,y: globals.dump(x, y+'T'), [ref,alt], [reffile, altfile])
	return complsnps

	
def checkRef(name):
	reffile = name +'RefT'
	altfile = name + 'AltT'
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = True)
	print "Loaded hg19"
	ref = globals.json(reffile)
	print "Loaded Ref"
	alt = globals.json(altfile)
	print "Loaded Alt"
	flip = []
	errors = []
	keys = ref.keys()
	for snppos in keys:
		t = snppos.split('pos')
		hg19snp = str(hg19[t[0]][int(t[1])-1]).upper()
		refsnp = ref[snppos].upper()
		altsnp = alt[snppos].upper()
		if hg19snp == refsnp:
			continue
		elif hg19snp == altsnp:
			ref[snppos] = hg19snp
			alt[snppos] = ref[snppos]
			flip.append(snppos)
		else:
			print "Error: Neither Ref nor Alt of SNP corresponds to hg19 sequence"
			errors.append(snppos)
	globals.dump(ref, reffile+'flipped')
	globals.dump(alt, altfile+'flipped')
	globals.dump(flip, name+'flips')
	globals.dump(errors, name+'errors')
	return [flip, errors]
		
def flipGeno(genofile, flip):
	lines = open(genofile).readlines()
	newgeno = open(genofile+'flipped', 'w')
	newgeno.write(lines[0])
	for l in lines[1:]:
		t = l.split('\t')
		if t[0] in flip:
			newg = map(lambda x: 2-int(x), t[1].split(','))
			newl = t[0] +'\t' 
			newl = reduce(lambda x,y: x+str(y) + ',', [newl]+newg)
			newgeno.write(newl.strip('\n'))
		else:
			newgeno.write(l)
			
#a = checkRef('../genotypes/CEUlowcov')
#flipGeno('../genotypes/CEUlowcovGeno', a[0])
names = ['../genotypes/CEUlowcov','../genotypes/YRIlowcov','../genotypes/CHBJPTlowcov','../genotypes/YRItrio', '../genotypes/CEUtrio']	
"""
for n in names:
	flips = globals.json(n+'flips')
	flipGeno(n, flips)	
"""

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
	globals.dump(ref, outputname+'Ref')
	globals.dump(alt, outputname+'Alt')
	
# Tests
def test():
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = True)
	for i in range(1,1000000):
		b = str(hg19['chr1'][i])
		
def testRef():
	Reffile = '../genotypes/CEUlowcovRef'
	Ref = simplejson.load(open(Reffile))
	keys = Ref.keys()
	for snppos in keys:
		continue