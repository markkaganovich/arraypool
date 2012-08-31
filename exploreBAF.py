import simplejson

report = 'Array1'

file = open(report)
arraylines = file.readlines()
file.close()

complement=  {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# correct the REF:ALT issue
# run through each snp on array, only keep the good ones, correct all the possible REF:ALT ones, 
refalt = simplejson.load(open('refalt'))
header = "SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,Allele1 - Plus,Allele2 - Plus,Chr,Position,SNP,Theta,R,X,Y,X Raw,Y Raw,B Allele Freq"
t = header.split(',')
snpi = t.index("SNP")
chri = t.index("Chr")
posi = t.index("Position")
Yi = t.index("Y")
Ri = t.index("R")

output = open(report.strip('.txt') + 'correctedYR','w')
for l in arraylines[11:]:
    t = l.split('\t')
    snp = t[snpi]
    a = snp.split('/')
    alleleA = a[0][1]
    alleleB = a[1][0]
    pos = 'chr'+t[chri]+'pos'+t[posi]
    print pos
    try:
		ra = refalt[pos]
		if ra['ref'] == alleleA and ra['alt'] == alleleB:  
		#	output.write(pos + '\t' +str(float(t[Yi])/float(t[Ri])) +'\n')
			output.write(str(float(t[Yi])/float(t[Ri])) +'\n')
		elif ra['ref'] == alleleB and ra['alt'] == alleleA:
			#output.write(pos + '\t' + str(1-float(t[Yi])/float(t[Ri]))+'\n' )
			output.write(str(1-float(t[Yi])/float(t[Ri]))+'\n' )	
    except KeyError:
        print "error"
