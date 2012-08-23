import simplejson

file = open('./corrected')
corrected = simplejson.load(file)
file.close()

file = open('./flipped')
flipped = simplejson.load(file)
file.close()

arrayBAFfile = open('./arrayBAFfile2','w')

arrayBAF = []
arraychosenIDs = []

omniidindexdic = simplejson.load(open('./Omni19IDposindexdic'))
refalt1indsfile = open('./refalt1indsfile','w')
refalt2indsfile = open('./refalt2indsfile','w')
relevantIDsfile = open('./relevantIDsfile2','w')
for c in corrected:
    t = c.split(',')
    #arrayBAF.append(float(t[5]))
    try:
    	refalt1indsfile.write(str(omniidindexdic[t[0]]+1)+'\n')
	arrayBAFfile.write(str(1-float(t[5])) +'\n')
	relevantIDsfile.write(str(t[0])+'\n')
    except KeyError:
	print t[0]

for f in flipped:
    t = f.split(',')
    #arrayBAF.append(1-float(t[5]))
    try:
    	refalt2indsfile.write(str(omniidindexdic[t[0]]+1)+'\n')
	arrayBAFfile.write(str(float(t[5])) +'\n')
	relevantIDsfile.write(str(t[0])+'\n')
    except KeyError:
	print t[0]
#file = open('./arrayBAF','w')
#simplejson.dump(arrayBAF, file)
#file.close()

#file = open('./arraychosenIDs','w')
#simplejson.dump(arraychosenIDs, file)
#file.close()

arrayBAFfile.close()
