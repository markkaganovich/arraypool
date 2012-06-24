exceptions =["chr4pos59480272", "chrXpos76056353", "chr15pos81350958", "chr7pos141666448", "chr2pos240488647", "chr4pos59481681", "chr4pos59482518"]

excinds = [205794, 347100, 750844, 783753, 819143, 879472,  1082264]

import simplejson

#file = open('omnisnpsHASH')
#hash = simplejson.load(file)
#file.close()
#snps = hash.keys()

#for s in range(0, len(snps)):
#    if snps[s] in exceptions:
#        print s


#need to open SnpMap instead of Human1M file
file = open('snpmap19')
snps19IDs = simplejson.load(file)
file.close()

snps18id = {}

file = open('omni19back')
omni19back = simplejson.load(file)
file.close()

file = open('omni18snps')
omni18 = simplejson.load(file)
file.close()

for s in range(0, len(omni19back)):
    if s in excinds:
        t = s+1
    else:
        t = s
    try:
        snps18id[snps19IDs[omni19back[s]]] = omni18[t]
    except KeyError:
        print omni19back[s]
    print s

file = open('snps18id','w')
simplejson.dump(snps18id, file)
file.close()
