import simplejson

newpool = './MK3.txt'

# import equimolar pool array inds and IDs

file = open('./relevantIDsfile')
lines = file.readlines()
file.close()

equalpoolIDs = []
for l in lines:
    equalpoolIDs.append(l.strip('\n'))


file = open(newpool)
lines = file.readlines()
file.close()

file = open(newpool + 'dic', 'w')

for l in lines[11:]:
    t = l.split(',')
    try:
	i = equalpoolIDs.index(t[0])
	file.write(str(i) + '\t' + l)
    except ValueError:
	#print "why"
	#print t
	continue

#file = open(newpool +'dic', 'w')
#simplejson.dump(newarraydic, file)
#file.close()
 

file = open('./omni19IDs')
omniIDs = simplejson.load(file)
file.close()

complement=  {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# correct the REF:ALT issue
# run through each snp on array, only keep the good ones, correct all the possible REF:ALT one

refalt = simplejson.load(open('refalt'))

correctedarraylines = []
notyet = []
flipped = []
correctedcomp = []
flipcompl = []
#dic = set(refalt.keys())
file = open(newpool + 'dic')
lines = file.readlines()
for l in lines:
    a= l.split('\t')
    b = a[1]
    t = b.split(',')
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
simplejson.dump(correctedarraylines, open(newpool +'corrected','w'))
simplejson.dump(flipped, open(newpool + 'flipped','w'))
simplejson.dump(correctedcomp, open(newpool + 'correctedcompl','w'))
simplejson.dump(flipcompl, open(newpool + 'flippedcoml','w'))
simplejson.dump(notyet, open(newpool + 'notyet','w'))
 
    

