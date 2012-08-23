import simplejson
'''
names = ['../genotypes/19CEUlowcov','../genotypes/19CHBJPTlowcov', '../genotypes/19YRIlowcov', '../genotypes/19CEUtrio', '../genotypes/19YRItrio']

for n in names:
	file = open(n+'SnpPos')
	snps = simplejson.load(file)
	file.close()

	snpdic = {}
	for x in range(0,len(snps)):
    		snpdic[snps[x]] =x

	file = open(n+'SnpPosDic','w')
	simplejson.dump(snpdic, file)
	file.close()

'''
#make reverse dic
'''
omniids = simplejson.load(open('./omni19IDs'))
omniposIDdic = {}
for k in omniids.keys():
    omniposIDdic[omniids[k]] = k

simplejson.dump(omniposIDdic, open('./omniposIDdic','w'))
'''

# consolidate RefAlt positions files
file = open('./omni19sorted')
omnisnps = simplejson.load(file)
file.close()

names =  ['../genotypes/19CEUlowcov','../genotypes/19CHBJPTlowcov', '../genotypes/19YRIlowcov', '../genotypes/19CEUtrio', '../genotypes/19YRItrio']

os = set(omnisnps)
refaltdic = {}

for name in names:
    file = open(name+'flatRefAlt')
    lines = file.readlines()
    for l in lines:
	t = l.split('\t')
	if t[0] in os:
	    refaltdic[t[0]] = {'ref': t[1], 'alt' : t[2].strip('\n')}
	
file = open('refalt','w')
simplejson.dump(refaltdic, file)
file.close()

