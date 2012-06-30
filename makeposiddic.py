import simplejson

omnisorted = simplejson.load(open('./omni19sorted'))
omniids = simplejson.load(open('./omni19IDs'))


idposindexdic = {}
#for k in omniids.keys():
 #   idposindexdic[k]= omnisorted.index(omniids[k])

for i in range(0, len(omnisorted)):
    snp = omnisorted[i]
    key = (key for key,value in omniids.items() if value==snp).next()
    idposindexdic[key] = i



file = open('./Omni19IDposindexdic','w')
simplejson.dump(idposindexdic, file)
file.close()

