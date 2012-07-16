import simplejson

reportfile = 'shortReport1.txt'


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
#file = open('omni19IDs')
#omni19IDs = simplejson.load(file)
#file.close()

file = open(reportfile)
lines = file.readlines()
file.close()

omni19results={}
del lines[0:11]
for l in lines:
    tokens = l.split('\t')
    if tokens.__len__() < 3:
        print l
        print "WHAT THE HELL"
    else:
        try:
	    if tokens[3] != 'NaN':
            	omni19results['chr'+tokens[1]+'pos'+tokens[2]] = float(tokens[3])
	    else:
            	omni19results['chr'+tokens[1]+'pos'+tokens[2]] = float(.00001)
		print tokens[0]
        except KeyError:
            print tokens[0]

file = open('omni19results'+reportfile, 'w')
simplejson.dump(omni19results, file)
file.close()

'''
omnisnps = simplejson.load(open('omni19sorted'))
nonzeroinds = simplejson.load(open('nonzeroinds'))

arrayfreq= []
file = open('arrayBfreq'+reportfile, 'w')
somni = set(nonzeroinds)
for i,s in enumerate(omnisnps):
    if i in somni:
	try: 
	    arrayfreq.append(omni19results[s])
	   # if i == 99404 or i == 99405:
	#	key = (key for key,value in omni19IDs.items() if value==s).next()
	#	print key
            file.write(str(omni19results[s])+'\n')
   	except KeyError:
            print "not in array"
            arrayfreq.append(0)
            file.write('0'+'\n')

file.close()
'''


