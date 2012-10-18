import simplejson

report = 'Array1'

arraysnps = map(lambda x: x.strip('\n'), open('omniexpresssnps').readlines())
arraysnpset = set(arraysnps)
file = open(report)
arraylines = file.readlines()
file.close()

complement=  {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# correct the REF:ALT issue
# run through each snp on array, only keep the good ones, correct all the possible REF:ALT ones, 
Ref = simplejson.load(open('Refomniexpress'))
Alt = simplejson.load(open('Altomniexpress'))
header = "SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,Allele1 - Plus,Allele2 - Plus,Chr,Position,SNP,Theta,R,X,Y,X Raw,Y Raw,B Allele Freq"
t = header.split(',')
snpi = t.index("SNP")
chri = t.index("Chr")
posi = t.index("Position")
Yi = t.index("Y")
Ri = t.index("R")

values = len(arraysnps) * [0]
output = open(report.strip('.txt') + 'correctedYR','w')
for l in arraylines[11:]:
	t = l.split('\t')
	snp = t[snpi]
	a = snp.split('/')
	alleleA = a[0][1]
	alleleB = a[1][0]
	pos = 'chr'+t[chri]+'pos'+t[posi]
	print pos
	if pos in arraysnpset:
		index = arraysnps.index(pos)
	else:
		continue
	try:
		r = Ref[pos]
		a = Alt[pos]
		if r == alleleA and a == alleleB:  
		#	output.write(pos + '\t' +str(float(t[Yi])/float(t[Ri])) +'\n')
			#output.write(str(float(t[Yi])/float(t[Ri])) +'\n')
			values[index] = float(t[Yi])/float(t[Ri])
		elif r == alleleB and a == alleleA:
			#output.write(pos + '\t' + str(1-float(t[Yi])/float(t[Ri]))+'\n' )
			#output.write(str(1-float(t[Yi])/float(t[Ri]))+'\n' )	
			values[index] = 1 - float(t[Yi])/float(t[Ri])
	except KeyError:
		print "error"

for v in values:
	output.write(str(v) + '\n')
