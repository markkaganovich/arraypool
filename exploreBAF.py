import simplejson

file = open('./MK1.txt')
arraylines = file.readlines()
file.close()


file = open('./omni19IDs')
omniIDs = simplejson.load(file)
file.close()

complement=  {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# correct the REF:ALT issue
# run through each snp on array, only keep the good ones, correct all the possible REF:ALT ones, 
# store good ones that aren't in teach Geno RefAlt file until the next one is opened

names = ['../genotypes/19CEUlowcov','../genotypes/19CHBJPTlowcov', '../genotypes/19YRIlowcov', '../genotypes/19CEUtrio', '../genotypes/19YRItrio']

correctedarraylines = []
notyet = []
flipped = []
correctedcomp = []
flipcompl = []
file = open(names[0]+'RefAlt')
refalt = simplejson.load(file)
file.close()
#dic = set(refalt.keys())
for l in arraylines[11:]:
    t = l.split(',')
    id = t[0]
    print id
    if t[9] != '-' and t[10] != '-':
	alleleA = t[9]
	alleleB = t[10]
	pos = omniIDs[id]
	print pos
	try:
	    ra = refalt[pos]
	    if ra['ref'] == alleleA and ra['alt'] == alleleB:
		    correctedarraylines.append(l)
	    elif ra['ref'] == alleleB and ra['alt'] == alleleA:
		    flipped.append(l)
	    elif ra['ref'] == complement[alleleB] and ra['alt'] == complement[alleleA]:
		    correctedcomp.append(l)
	    elif ra['ref'] == complement[alleleA] and ra['alt'] == complement[alleleB]:
		    flipcompl.append(l)
	    else:
		    notyet.append(l)
	except KeyError:
	    print "error"
'''
simplejson.dump(correctedarraylines, open('corrected','w'))
simplejson.dump(flipped, open('flipped','w'))
simplejson.dump(correctedcomp, open('correctedcompl','w'))
simplejson.dump(flipcompl, open('flippedcoml','w'))
simplejson.dump(notyet, open('notyet','w'))
'''
for name in names[1:]:
    if notyet == []:
	break
    else:
	refalt = simplejson.load(open(name+'RefAlt'))
        dic = set(refalt.keys())
    for ny in notyet:
	t = ny.split(',')
	id = t[0]
	print id
	alleleA = t[9]
	alleleB = t[10]
	pos = omniIDs[id]
	if pos in dic:
	    if refalt[pos]['ref'] == alleleA and refalt[pos]['alt'] == alleleB:
		    correctedarraylines.append(l)
	    elif refalt[pos]['ref'] == alleleB and refalt[pos]['alt'] == alleleA:
		    flipped.append(l)
	    elif refalt[pos]['ref'] == complement[alleleB] and refalt[pos]['alt'] == complement[alleleA]:
		    correctedcompl.append(l)
	    elif refalt[pos]['ref'] == complement[alleleA] and refalt[pos]['alt'] == complement[alleleB]:
		    flipcompl.append(l)
	    else:
		newny.append(l)
	
    notyet = newny
	    
simplejson.dump(correctedarraylines, open('corrected','w'))
simplejson.dump(flipped, open('flipped','w'))
simplejson.dump(correctedcomp, open('correctedcompl','w'))
simplejson.dump(flipcompl, open('flippedcoml','w'))
simplejson.dump(notyet, open('notyet','w'))

