import simplejson
from pygr import worldbase


def checkRef(Reffile, Altfile):
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = True)
	Ref = simplejson.load(open(Reffile))
	Alt = simplejson.load(open(Altfile))
	flip = []
	errors = []
	keys = Ref.keys()
	for snppos in keys:
		t = snppos.split('pos')
		hg19snp = str(hg19[t[0]][int(t[1])-1]).upper()
		refsnp = Ref[snppos].upper()
		altsnp = Alt[snppos].upper()
		if hg19snp == refsnp:
			continue
		elif hg19snp == altsnp:
			Ref[snppos] = hg19snp
			Alt[snppos] = Ref[snppos]
			flip.append(snppos)
		else:
			print "Error: Neither Ref nor Alt of SNP corresponds to hg19 sequence"
			errors.append(snppos)
	file = open(Reffile+'flipped','w')
	simplejson.dump(Ref, file)
	file.close()
	file = open(Altfile+'flipped','w')
	simplejson.dump(Alt, file)
	file.close()
	return [flip, errors]
	
def flipGeno(genofile, flip):
	lines = open(genofile).readlines()
	newgeno = open(genofile+'flipped', 'w')
	newgeno.write(lines[0])
	for l in lines[1:]:
		t = l.split('\t')
		if t[0] in flip:
			newg = map(lambda x: 2-int(x), l.split(','))
			newl = t[0] +'\t' 
			map(lambda x: newl + str(x) + ',')
			newgeno.write(newl.strip('\n'))
		else:
			newgono.write(l)
			
a = checkRef('../genotypes/CEUlowcovRef', '../genotypes/CEUlowcovAlt')
flipGeno('../genotypes/CEUlowcovGeno', a[0])

#modified to the 19 version
def parse1KGvcf(vcffile, outputname):
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
	simplejson.dump(ref, outputref)
	simplejson.dump(alt, outputalt)