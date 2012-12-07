import simplejson
import os
import re

def jsonload(filename):
	file = open(filename)
	result = simplejson.load(file)
	file.close()
	return result

def jsondump(data, filename, ds='./'):
	n = len(filter(lambda x: re.findall(r"^"+filename+"[0-9]*",x), os.listdir(ds)))

	"""
	renaming existing file in an odd way, for development purposes; usually people make new file have
	the higher number, here its the old one
	"""
	if n > 0:
		os.rename(ds +filename, ds +filename+ str(n+1))	  
	file = open(filename,'w')
	simplejson.dump(data, file)
	file.close()

compl = {'G' :'C', 'C' : 'G', 'A' : 'T', 'T':'A'}
chrset = range(1,23) + ['X','Y']

class Ref():
	
	def __init__(self, name):
		self.current = name
		
	syn = {}
	syn['hg19'] = ['GRCh37']
	syn['GRCh37'] = ['hg19']
	syn['GRCh36'] = ['hg18']

