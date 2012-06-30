import simplejson

file = open('./corrected')
corrected = simplejson.load(file)
file.close()

file = open('./flipped')
flipped = simplejson.load(file)
file.close()

arrayBAF = []
arraychosenIDs = []

for c in corrected:
    t = c.split(',')
    arrayBAF.append(float(t[5]))
    arraychosenIDs.append(t[0])

for f in flipped:
    t = f.split(',')
    arrayBAF.append(1-float(t[5]))
    arraychosenIDs.append(t[0])

file = open('./arrayBAF','w')
simplejson.dump(arrayBAF, file)
file.close()

file = open('./arraychosenIDs','w')
simplejson.dump(arraychosenIDs, file)
file.close()

omniidindexdic = simplejson.load(open('./omni19iddic'))
relevantinds = []
for a in arraychosenIDs:
    relevantinds.append(omniidindexdic[a] + 1)

simplejson.dump(relevantinds, open('./relevantinds','w'))


