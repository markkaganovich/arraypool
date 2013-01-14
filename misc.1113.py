'''
misc code for arrayvcf.py that wasn't really used

'''
def in1kg(poollines, names):
	lines = []
	for n in names:
		lines.extend(getlines(n))
	result = []
	for l in poollines:
		if l in lines:
			result.append(True)
		else:
			result.append(False)
	return result

def getlines(name):
	vfile = open(name, 'r')
	vcf_reader = vcf.Reader(vfile)
	lines = vcf_reader.samples
	return lines

def comparearrays(array1, array2):
	thetai = 10
	snpidi = 0
	a1 = open(array1)
	a2 = open(array2)
	out = open(array1+array2 +'.diff.theta.R', 'w')
	for l1, l2 in zip(a1, a2):
			t1 = l1.split('\t')
			t2 = l2.split('\t')
			if t1[snpidi] == t2[snpidi] and t1[thetai] != 0:
				out.write(str((float(t1[thetai]) - float(t2[thetai]))/sum(float(t1[thetai]), float(t2[thetai]))/2) +'\n')
	out.close()
