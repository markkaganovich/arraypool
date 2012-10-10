import simplejson

def parsehapmapchrom(chrom, rshash, refalt={}):
	people = ['NA19140','NA19154','NA19173','NA19203','NA19206','NA19211','NA19222']
	filename = '../genotypes/hapmapchr'+str(chrom)
	hf = open(filename)
	lines = hf.readlines()
	header = lines[0].split(' ')
	out = open(filename+'genotype','w')
	posout = open(filename+'SnpPos','w')
	rsidi = header.index('rs#')
	chromi = header.index('chrom')
	posi = header.index('pos')	
	snpi = header.index('alleles')
	peoplewithgenos = filter(lambda x: x in header, people)
	print peoplewithgenos
	peoplei = map(lambda x: header.index(x), peoplewithgenos)
	refalt = {}
	refaltout = open(filename+'RefAlt', 'w')
	for l in lines[1:]:
		t = l.split(' ')
		rsid = t[rsidi]
		chrom = t[chromi]
		pos = t[posi]
		snps = t[snpi].split('/')
		ref = snps[0]
		alt = snps[1]
		for p in peoplei:
			t[p]
			genotypes = map(lambda x: t[x], peoplei)
			genocount = map(lambda x: len(filter(lambda y: alt == y, x)), genotypes)
			try:
				SnpPos = rshash[rsid]
				map(lambda x: out.write(str(x)+','), genocount[:-1])
				out.write(str(genocount[-1])+'\n')
				posout.write(SnpPos+'\n')
				refalt[SnpPos]={}
				refalt[SnpPos]['ref'] = ref
				refalt[SnpPos]['alt'] = alt
			except KeyError:
				pass
	simplejson.dump(refalt, refaltout)
	return refalt
			
file = open('./hapmap_rsid_hash_lines')
rsid2pos = {}
rlines = file.readlines()
for r in rlines:
	t = r.split('\t')
	rsid2pos[t[0]] = t[1].strip('\n')
file.close()
refalt = {}	
for c in range(1,23):
	refalt = parsehapmapchrom(c, rsid2pos, refalt)
file = open('../genotypes/hapmapRefAlt', 'w')
simplejson.dump(refalt, file)
file.close()
