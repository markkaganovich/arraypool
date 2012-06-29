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
for l in arraylines[11:]:
    t = l.split(',')
    id = t[0]
    if t[9] != '-' and t[10] != '-':
	    alleleA = t[9]
	    alleleB = t[10]
    if refalt[omniIDs[id]]['ref'] == alleleA and refalt[omniIDs[id]]['alt'] == alleleB:
	    correctedarraylines.append(l)
    elif refalt[omniIDs[id]]['ref'] == alleleB and refalt[omniIDs[id]]['alt'] == alleleA:
	    flipped.append(l)
    elif refalt[omniIDs[id]]['ref'] == complement[alleleB] and refalt[omniIDs[id]]['alt'] == complement[alleleA]:
	    correctedcompl.append(l)
    elif refalt[omniIDs[id]]['ref'] == complement[alleleA] and refalt[omniIDs[id]]['alt'] == complement[alleleB]:
	    flipcompl.append(l)
    else:
	    notyet.append(l)

simplejson.dump(correctedarraylines, open('corrected','w'))
simplejson.dump(flipped, open('flipped','w'))
simplejson.dump(correctedcompl, open('correctedcompl','w'))
simplejson.dump(flipcompl, open('flippedcoml','w'))
smplejson.dump(notyet, open('notyet','w'))


