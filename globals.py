import simplejson
import os

datasource = '/Users/markkaganovich/Dropbox/data/'

def json(filename, datasource = datasource):
    file = open(datasource + filename)
    result = simplejson.load(file)
    file.close()
    return result

def dump(data, filename, datasource = ''):
    if filename in os.listdir(datasource): 
        os.rename(datasource +filename, datasource +'temp/' +filename+'2')
    file = open(datasource + filename,'w')
    simplejson.dump(data, file)
    file.close()

compl = {'G' :'C', 'C' : 'G', 'A' : 'T', 'T':'A'}


'''
def addfeatures(classname, feature, hash, fun, args = None, data = None):
    methods = getattr(classname, hash)
    methods[feature] = {}
    methods[feature]['fun'] = fun
    if args != None:
        methods[feature]['args'] = args
    if data != None:
        methods[feature]['data'] = data
    dump(methods, classname.__name__+str(hash), './dbase/')
'''
'''
    if 'args' in kwargs.keys():
        methods['args'] = kwargs['args']
    if 'data' in kwargs.keys():
        methods['data'] = kwargs['data']
'''
'''
'''


