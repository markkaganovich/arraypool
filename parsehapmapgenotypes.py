import os

def parsehapmapchrom(chrom):
	""" inputs: rsid hash (hash of snp positions and rsIDs)
				hapmapchrN downloaded from hapmap3 site
				people: list of cell line IDs that we are looking for
	"""
	
	hapmapfile = '../genotypes/hapmapchr'+str(chrom)
	if 'rsid2poshash' not in os.listdir('./'):
		makerhash()
	else:
		rshash = glob.json('rsid2poshash')	
	people = ['NA19140','NA19154','NA19173','NA19203','NA19206','NA19211','NA19222']
		
	with open(hapmapfile) as f:
		lines = f.readlines()	
	header = lines[0].split(' ')
	
	rsidi = header.index('rs#')
	snpi = header.index('alleles')
	
	peoplewithgenos = filter(lambda x: x in header, people)
	print peoplewithgenos
	peoplei = map(lambda x: header.index(x), peoplewithgenos)
	
	refhash = {}
	althash = {}
	out = open(hapmapfile+'genotype','w')
	
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

def makerhash():
	rsid2pos = {}
	hapmaprsidfile = 'hapmap_rsid_hash_lines'
	assert hapmaprsidfile in os.listdir('./'), "Need this file {0}".format(hapmaprsidfile)			
	file = open('./hapmap_rsid_hash_lines')
	with open('./hapmap_rsid_hash_lines') as f:
		rlines = f.readlines()
	for r in rlines:
		t = r.split('\t')
		rsid2pos[t[0]] = t[1].strip('\n')	
	glob.dump(rsid2pos, 'rsid2poshash')

	

