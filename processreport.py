import os

#process report
def splitreport(f, dir):
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
				outputfile = '25M' + l.split('\t')[i]
			except IndexError:
				continue
			if l.split('\t')[i] != '3.17':
				continue
			else:
				if outputfile in os.listdir(dir):
					out.write(l)
				else:
					out = open(outputfile, 'w')
					out.write(l)
					print out