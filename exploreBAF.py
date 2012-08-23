import simplejson

report = 'MKReportbySNP2.txt'

file = open(report)
arraylines = file.readlines()
file.close()

complement=  {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# correct the REF:ALT issue
# run through each snp on array, only keep the good ones, correct all the possible REF:ALT ones, 
refalt = simplejson.load(open('refalt'))

output = open(report.strip('.txt') + 'correctedYR','w')
for l in arraylines[11:]:
    t = l.split('\t')
    id = t[0]
    print id
    #if t[9] != '-' and t[10] != '-':
    snp = t[14]
    a = snp.split('/')
    alleleA = a[0][1]
    alleleB = a[1][0]
    pos = 'chr'+t[12]+'pos'+t[13]
    print pos
    try:
	ra = refalt[pos]
	print ra
	print t[18]
        if ra['ref'] == alleleA and ra['alt'] == alleleB and t[18]!='-' and t[18]!='+':  
	    output.write('chr'+t[12]+'pos'+t[13] + '\t' +str(float(t[21])/float(t[19])) +'\n')
#t[18]+'\n')
	elif ra['ref'] == alleleB and ra['alt'] == alleleA and t[18]!='-' and t[18]!='+':
	    output.write('chr'+t[12]+'pos'+t[13]+'\t'+ str(1-float(t[21])/float(t[19]))+'\n' )
#str(1-float(t[18])) +'\n')
	#elif ra['ref'] == complement[alleleB] and ra['alt'] == complement[alleleA]:
	 #   correctedcomp.append(l)
	#elif ra['ref'] == complement[alleleA] and ra['alt'] == complement[alleleB]:
	 #   flipcompl.append(l)
	#else:
	 #   notyet.append(l)	
    except KeyError:
        print "error"
'''
simplejson.dump(correctedarraylines, open('corrected','w'))
simplejson.dump(flipped, open('flipped','w'))
simplejson.dump(correctedcomp, open('correctedcompl','w'))
simplejson.dump(flipcompl, open('flippedcoml','w'))
simplejson.dump(notyet, open('notyet','w'))
'''
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
'''
