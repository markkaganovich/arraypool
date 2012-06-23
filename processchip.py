import simplejson

reportfile = 'MK3.txt'

file = open(reportfile)
lines = file.readlines()
file.close()

del lines[0:11]

file = open('snps18id')
snps18id = simplejson.load(file)
file.close()

omni18results={}

for l in lines:
    tokens = l.split(',')
    if tokens.__len__() < 3:
        print l
        print "WHAT THE HELL"
    else:
        try:
            omni18results[snps18id[tokens[0]]] = float(tokens[5])
        except KeyError:
            print tokens[0]

file = open('omni18results'+reportfile, 'w')
simplejson.dump(omni18results, file)
file.close()


omnisnps = simplejson.load(open('omnisnpmap18'))
arrayfreq= []
file = open('arrayBfreq'+reportfile, 'w')
for s in omnisnps:
    try: 
        arrayfreq.append(omni18results[s])
        file.write(str(omni18results[s])+',')
    except KeyError:
        print "not in array"
        arrayfreq.append(0)
        file.write('0,')

file.close()



