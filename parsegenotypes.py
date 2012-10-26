import simplejson
from pygr import worldbase
import hashes

def filterSNPs(name):
	Reffile = name +'Ref'
	Altfile = name + 'Alt'
	Ref = simplejson.load(open(Reffile))
	print "Loaded Ref"
	print len(Ref)
	Alt = simplejson.load(open(Altfile))
	print "Loaded Alt"
	print len(Alt)
	keys = Ref.keys()
	complsnps = []
	for snppos in keys:
		if hashes.compl[Ref[snppos].upper()] == Alt[snppos].upper():
			complsnps.append(snppos)
			del Ref[snppos]
			del Alt[snppos]
	print len(Ref)
	print len(Alt)
	file = open(Reffile,'w')
	simplejson.dump(Ref, file)
	file.close()
	file = open(Altfile,'w')
	simplejson.dump(Alt, file)
	file.close()
	return complsnps

	
def checkRef(name):
	Reffile = name +'Ref'
	Altfile = name + 'Alt'
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = True)
	print "Loaded hg19"
	Ref = simplejson.load(open(Reffile))
	print "Loaded Ref"
	Alt = simplejson.load(open(Altfile))
	print "Loaded Alt"
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
	file = open(name+'flips','w')
	simplejson.dump(flip, file)
	file.close()
	file = open(name+'errors', 'w')
	simplejson.dump(errors, file)
	file.close()
	return [flip, errors]
	
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