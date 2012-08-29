import simplejson

def parsehapmapchrom(chrom):
	people = ['NA19140','NA19154', 'NA19145','NA19173'.'NA19203','NA19206','NA19211','NA19222']
	filename = '../genotypes/hapmapchr'+str(chrom)
	hf = open(filename)
	lines = hf.readlines()
	header = lines[0].split(' ')
	out = open(filename+'genotype','w')
	rsidout = open(filename+'posrsid','w')
	rsidi = header.index('rs#')
	chromi = header.index('chrom')
	posi = header.index('pos')	
	snpi = header.index('alleles')
	peoplei = map(lambda x: header.index(x), people)
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
			print genocount
			map(lambda x: out.write(str(x)+','), genocount[:-1])
			out.write(str(genocount[-1])+'\n')
			rsidout.write(rsid)
	
for c in range(1,23):
	parsehapmapchrom(c)