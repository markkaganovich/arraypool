import simplejson

file = open('./omni19resultsshortReport1.txt')
results1 = simplejson.load(file)
file.close()

file = open('./omni19resultsshortReport3.txt')
results2 = simplejson.load(file)
file.close()
'''
set1 = set(results1.keys())
set2 = set(results2.keys())

if set1 == set2:
    print "same"
else:
    print "different"

file = open('./omni19snps','w')
simplejson.dump(results1.keys(), file)
file.close()
'''

output1 = open('./arrayequaltheta','w')
output2 = open('./arrayspikedtheta','w')
for k in results1.keys():
    output1.write(str(results1[k]) + '\n')
    output2.write(str(results2[k]) + '\n')

output1.close()
output2.close()





