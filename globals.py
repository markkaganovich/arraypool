import simplejson
import os

datasource = '/Users/markkaganovich/Dropbox/data/'

def json(filename, datasource = ''):
    file = open(datasource + filename)
    result = simplejson.load(file)
    file.close()
    return result

def dump(data, filename, datasource = './'):
    #if filename in os.listdir(datasource): 
    #    os.rename(datasource +filename, datasource +'temp/' +filename+'2')
    file = open(filename,'w')
    simplejson.dump(data, file)
    file.close()

compl = {'G' :'C', 'C' : 'G', 'A' : 'T', 'T':'A'}



