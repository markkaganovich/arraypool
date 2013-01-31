import os

#process report
def splitreport(f, arrayname):
	dir = './'
	r = open(f)
	header = None
	for l in r:
		if "[Data]" in l:
			header = r.next().split('\t')
			print header
			continue
		if header:
			i = header.index('Sample ID')
			try:
				outputfile = arrayname + l.split('\t')[i]
			except IndexError:
				continue
			if outputfile in os.listdir(dir):
				out.write(l)
			else:
				out = open(outputfile, 'w')					
				out.write(l)
				print out