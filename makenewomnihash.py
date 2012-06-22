import simplejson

'''
file = open('./SNP_Map.csv')
lines = file.readlines()
file.close()
file = open('snpmapliftoverinput','w')

del lines[0]
for l in lines:
    t = l.split(',')
    file.write('chr'+t[2] +  ':' + t[3]+'-'+str(int(t[3])+1)+'\n')

file.close()
'''
lines = open('./hglft_snpmap.bed.txt').readlines()
omnisnpmap18 = []
for l in lines:
    t = l.split(':')
    omnisnpmap18.append(t[0]+'pos'+t[1].split('-')[0])

