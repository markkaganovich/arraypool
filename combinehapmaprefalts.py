import simplejson

#combine hapmapchr(i)RefAlt into one dictionary hapmapRefAlt
def combine():
	hapmapRefAlt = {}
	for c in range(1,23):
		f = open('./hapmapchr'+str(c)+'RefAlt')
		hmrefalt = simplejson.load(f)
		f.close()
		hapmapRefAlt.update(hmrefalt)
		outf = open('./hapmapRefAlt','w')
		simplejson.dump(hapmapRefAlt, outf)	
	

names = ['./19CEUlowcov','../genotypes/19CHBJPTlowcov', '../genotypes/19YRIlowcov', '../genotypes/19CEUtrio', '../genotypes/19YRItrio']
	
f = open('./hapmapRefAlt')
hapmapRefAlt = simplejson.load(f)
f.close()

f = open(names[0]+'RefAlt')
ra = simplejson.load(f)
f.close()

for k in hapmapRefAlt.keys():
	try:
		if ra[k] == hapmapRefAlt[k]:
			matches.append(k)
	except KeyError:
		pass
		

