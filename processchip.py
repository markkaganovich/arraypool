import simplejson

reportfile = 'MK3.txt'


#file = open('./omni19sorted')
#omnisnps = simplejson.load(file)
#file.close()
'''
file = open('./SNP_Map.csv')
lines = file.readlines()
file.close()
del lines[0]
omni19IDs ={}
for l in lines:
    t = l.split(',')
    omni19IDs[t[1]] = 'chr'+t[2]+'pos'+t[3]
'''
file = open('omni19IDs')
omni19IDs = simplejson.load(file)
file.close()

file = open(reportfile)
lines = file.readlines()
file.close()

omni19results={}
del lines[0:11]
for l in lines:
    tokens = l.split(',')
    if tokens.__len__() < 3:
        print l
        print "WHAT THE HELL"
    else:
        try:
            omni19results[omni19IDs[tokens[0]]] = float(tokens[5])
        except KeyError:
            print tokens[0]

file = open('omni19results'+reportfile, 'w')
simplejson.dump(omni19results, file)
file.close()


omnisnps = simplejson.load(open('omni19sorted'))


arrayfreq= []
file = open('arrayBfreq'+reportfile, 'w')
for s in omnisnps:
    try: 
        arrayfreq.append(omni19results[s])
        file.write(str(omni19results[s])+'\n')
    except KeyError:
        print "not in array"
        arrayfreq.append(0)
        file.write('0'+'\n')

file.close()



