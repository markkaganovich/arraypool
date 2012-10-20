import simplejson
from pygr import worldbase


def checkRef(Ref, Alt):
	hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19(download = TRUE)
	Refcorrected = {}
	Altcorrected = {}
	flip = []
	for snppos in Ref.keys():
		t = snppos.split('pos')
		hg19snp = str(hg19[t[0]][int(t[1])-1])
		if hg19snp == Ref[snppos]:
			Refcorrected[snppos] = Ref[snppos]
			Altcorrected[snppos] = Alt[snppos]
		elif hg19snp = Alt[snppos]:
			Refcorrected[snppos] = Alt[snppos]
			Altcorrected[snppos] = Ref[snppos]
			flip.append(snppos)
		else:
			print "Error: Neither Ref nor Alt of SNP corresponds to hg19 sequence"
	return [Refcorrected, Altcorrected, flip]
			
a = checkRef('../genotypes/CEUlowcovRef', '../genotypes/CEUlowcovAlt')


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