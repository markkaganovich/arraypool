import simplejson

reportfile = 'MK1.txt'

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
    try:
        omni18results[snps18id[tokens[0]]] = float(tokens[5])
    except KeyError:
        print tokens[0]

file = open('omni18results', 'w')
simplejson.dump(omni18results, file)
file.close()

