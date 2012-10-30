import simplejson

def parsehapmapchrom(chrom, rshash):
	people = ['NA19140','NA19154','NA19173','NA19203','NA19206','NA19211','NA19222']
	filename = '../genotypes/hapmapchr'+str(chrom)
	with open(filename) as f:
		lines = f.readlines()
	
	header = lines[0].split(' ')
	
	out = open(filename+'genotype','w')
	
	rsidi = header.index('rs#')
	#chromi = header.index('chrom')
	#posi = header.index('pos')	
	snpi = header.index('alleles')
	
	peoplewithgenos = filter(lambda x: x in header, people)
	print peoplewithgenos
	peoplei = map(lambda x: header.index(x), peoplewithgenos)
	
	refhash = {}
	althash = {}
	
	for l in lines[1:]:
		t = l.split(' ')
		rsid = t[rsidi]
		snps = t[snpi].split('/')
		ref = snps[0]
		alt = snps[1]
		
		for p in peoplei:
			genotypes = map(lambda x: t[x], peoplei)
			genocount = map(lambda x: len(filter(lambda y: alt == y, x)), genotypes)
			try:
				snppos = rshash[rsid]
				newline = snppos + '\t'
				newline = reduce(lambda x,y: x+str(y) + ',', [newline] + genocount)
				out.write(newline+'\n')

				refhash[snppos] = ref
				althash[snppos] = alt
			except KeyError:
				pass
	out.close()
	return [refhash, althash]

rsid2pos = {}			
file = open('./hapmap_rsid_hash_lines')
with open('./hapmap_rsid_hash_lines') as f:
	rlines = f.readlines()

for r in rlines:
	t = r.split('\t')
	rsid2pos[t[0]] = t[1].strip('\n')

ref = {}
alt = {}	
genotype = open('hapmapGeno','w')
for c in range(1,23):
	[r,a] = parsehapmapchrom(c, rsid2pos)
	ref.update(r)
	alt.update(a)
	with open('../genotypes/hapmapchr'+str(c)) as g:
		lines = g.readlines()
		map(lambda l: genotype.write(l), lines)
genotype.close()
globals.dump(ref, '../genotypes/hapmapRef')
globals.dump(alt, '../genotypes/hapmapAlt')

