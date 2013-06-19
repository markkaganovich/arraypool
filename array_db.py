#array_db.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table, select
import sys
import json
import commands
import os
import csv 
import re
import numpy as np


db = create_engine('sqlite:///../cancergenomes/GENOTYPES_v2.db', echo = False)
db_hapmap = create_engine('sqlite:///../cancergenomes/GENOTYPES.db')
db_array = create_engine('sqlite:///../arraydata/arrays.db')


Session = sessionmaker(db)
session = Session()
metadata = MetaData(db)

Session_array = sessionmaker(db_array)
session_array = Session_array()
metadata_array = MetaData(db_array)

metadata_hapmap = MetaData(db_hapmap)

hapmap_table_name = "hapmap"
array_table_name = "Array1S_finalreport"

def copy_my_table(name, table):
    args = []
    for c in table.columns:
        args.append(c.copy())
    for c in table.constraints:
        args.append(c.copy())
    return Table(name, table.metadata, *args)

def get_rsids(name, table):
    if name in os.listdir('./'):
        print "loading from file"
        rsids = json.load(open(name, 'r'))
    else:
        rsids = set([])
        s = table.select()
        rs = s.execute()
        for row in rs:
            #rsids.add(getattr(row, 'snp name'))
            rsids.add(getattr(row, 'rsid'))     
        json.dump(list(rsids), open(name, 'w'))
    return rsids

def get_samples(table, gtype = 'kg'):
    if gtype == 'kg':
        samples = table.columns.keys()
    if gtype == 'hapmap':
        samples = table.columns.keys()
    return samples


#kg_table = Table(kg_table_name, metadata, autoload = True)
hapmap_table = Table(hapmap_table_name, metadata_hapmap, autoload = True)
array_table = Table(array_table_name, metadata_array, autoload = True)
print "post-table"



kg_files = ['CEU.low_coverage.2010_09.genotypes.vcf', 'YRI.low_coverage.2010_09.genotypes.vcf', 'CHBJPT.low_coverage.2010_09.genotypes.vcf', 'CEU.trio.2010_09.genotypes.vcf', 'YRI.trio.2010_09.genotypes.vcf']

class KG:
    def __init__(self, name, filename, trio_bool):
        self.name = name
        self.lines = open('../1000GenomesData/'+filename).readlines()
        if trio_bool:
            self.header = self.lines[14].split('\t')
            self.lines = self.lines[15:]
        else: 
            self.header = self.lines[17].split('\t')
            self.lines = self.lines[18:]
        self.rs_index = json.load(open(name+'_rs_index'))

print "loading files...."
ceu = KG('ceu', kg_files[0], 0)
print "ceu"
yri = KG('yri', kg_files[1], 0)
chbjpt = KG('chbjpt', kg_files[2], 0)
ceu_trio = KG('ceu_trio', kg_files[3], 1)
yri_trio = KG('yri_trio', kg_files[4], 1)

kgs = [ceu, yri, chbjpt, ceu_trio, yri_trio]

kg_rsids = []
for kg in kgs:
    for l in kg.lines:
        rs = l.split('\t')[kg.header.index('ID')]
        kg_rsids.append(rs)

kg_rsids = set(kg_rsids)
hapmap_rsids = get_rsids('hapmap_rsids', hapmap_table)
array_rsids = get_rsids('array_rsids', array_table)

inboth = kg_rsids.intersection(hapmap_rsids)
inall = list(inboth.intersection(set(array_rsids)))
select_rsids = set(inall[1:1000])

def make_rsindex(kg_object):
    rs_index = {}
    for i,l in enumerate(kg_object.lines):
        rs_index[l.split('\t')[kg_object.header.index('ID')]] = i
    json.dump(rs_index, open(kg_object.name+'_rs_index', 'w'))





def get_rows(attr = 'rs#', selected_attr = [], table = None):
    rows = []
    s = table.select()
    rs = s.execute()
    for r in rs:
        if getattr(r, attr) in selected_attr:
            rows.append(r)
    return rows

#kg_rsid_rows = get_rows(attr = 'rs#', selected_attr = select_rsids, table = kg_table)


'''
for i in range(1,len(hapmap_rsid_rows)-1): 
    if getattr(hapmap_rsids_sorted[i-1], 'rs#') == getattr(hapmap_rsids_sorted[i], 'rs#') \
        or (getattr(hapmap_rsids_sorted[i+1], 'rs#') == getattr(hapmap_rsids_sorted[i], 'rs#')):
        hapmap_rsid_rows.remove(hapmap_rsids_sorted[i])
'''
'''
given list of cell line IDs, find if they are coming from hapmap or 1kg, 
get their GENOTYPES for all the rsids 
'''
pool1_samples = ["NA18516", "NA18517", "NA18579", "NA18592", "NA18561", "NA07357", "NA06994", "NA18526", "NA12004", "NA19141", "NA19143", "NA19147", "NA19152", "NA19153", "NA19159", "NA19171", "NA19172", "NA19190", "NA19207", "NA19209", "NA19210", "NA19225", "NA18856", "NA18858", "NA18562", "NA18563", "NA18853", "NA18861"]


hapmap_samples = get_samples(table = hapmap_table, gtype = 'hapmap')


pool2 = json.load(open('pool2'))
pool_samples = pool2


def find_rs_line_kg(k_object, rs, s):
    geno = {}
    geno['source'] = k_object.name
    try:
        index = k_object.rs_index[rs]
    except KeyError:
        geno['genotype'] = 0
        return geno
    l = k_object.lines[index]
    if rs == l.split('\t')[k_object.header.index('ID')]:  
        info = l.split('\t')[k_object.header.index(s)]
        g = info.split(':')[0]
        g_split = re.split('[/|\\\.|]', g)
        if len(g_split) == 2:
            genotype = sum(map(lambda x: int(x), g_split))
            geno['ref'] = l.split('\t')[k_object.header.index('REF')]        
            geno['alt'] = l.split('\t')[k_object.header.index('ALT')]
            geno['genotype'] = genotype
            return geno

def find_rs_line_hapmap(rs, s):
    geno = {}
    for h in hapmap_rsid_rows:
            if getattr(h, 'rs#') == rs:
                geno = {}
                geno['source'] = 'hapmap'
                geno['ref'] = h[2].split('/')[0]
                geno['alt'] = h[2].split('/')[1]
                s_ind = hapmap_samples.index(s.lower())
                g = h[s_ind]
                geno['genotype'] = 2-len(filter(lambda x: geno['ref'] == x, g))
                print geno
                return geno
    geno['genotype'] = 0
    return geno



#kg_rsids_sorted = sorted(kg_rsid_rows, key=lambda r: getattr(r, 'rs#'))
#hapmap_rsids_rows_sorted = sorted(hapmap_rsid_rows, key=lambda r: getattr(r, 'rs#'))

rs = list(inall)[0:20000]
rs_lines = []
s = array_table.select()
r = s.execute()
setrs = set(rs)
hapmap_rsid_rows = get_rows(attr = 'rs#', selected_attr = setrs, table = hapmap_table)
rows = filter(lambda x: getattr(x, 'snp name') in setrs, r)

def find_rs_line_array(rs, rows):
    for r in rows:
        if getattr(r, 'snp name') == rs:
            print r
            return getattr(r, 'b allele freq')



# will need to split this up by specific array sample



g_rs = {}
for r in rs:
    genotypes = {}
    for p in pool_samples:
        #print p
        for k in kgs:        
            if p in k.header:
                if p == pool_samples[3]:
                    print k.name
                    print "here"
                genotypes[p] = find_rs_line_kg(k, r, p)        
        if p.lower() in hapmap_samples:
            if p == pool_samples[3]:
                    print "hapmap"
            genotypes[p] = find_rs_line_hapmap(r, p)
            print genotypes[p]
    g_rs[r]=genotypes


needtoflip = []
hapmaprefs = {}
kgrefs = {}
delete = []
for r in rs:
    hapmaprefs[r] = []
    kgrefs[r] = []
    for k in g_rs[r]:
        if g_rs[r][k]['source'] == 'hapmap':
            try:
                hapmaprefs[r].append(g_rs[r][k]['ref'])
            except KeyError:
                continue
        else:
            try:
                kgrefs[r].append(g_rs[r][k]['ref'])
            except KeyError:
                continue
    if len(kgrefs[r]) < 1:
        delete.append(r)
        continue
    if len(set(hapmaprefs[r])) != 1:
        print "error in hapmap refs"
    else:
        if hapmaprefs[r][0] != kgrefs[r][0]:
            needtoflip.append(r)

g_rs = {key: value for key,value in g_rs.items() if key not in delete}

for r in needtoflip:
    for k in g_rs[r]:
        if g_rs[r][k]['source'] == 'hapmap':
            tempref = g_rs[r][k]['ref']
            g_rs[r][k]['ref'] = g_rs[r][k]['alt']
            g_rs[r][k]['alt'] = tempref
            g_rs[r][k]['genotype'] = 2- int(g_rs[r][k]['genotype'])


        #tempref = hapmaprefs[r].values()[0]

    




g_freqs = []
g = np.empty(len(g_rs.keys())* len(pool_samples)).reshape(len(g_rs.keys()), len(pool_samples))
for r,i in enumerate(g_rs.keys()):
    #g_freqs.append(sum(map(lambda x: g_rs[r][x]['genotype'], g_rs[r].keys()))/ float(len(g_rs[r].keys()*2)))
    for s,j in enumerate(pool_samples):
        g[r,s] = g_rs[i][j]['genotype']



a_freqs = []
for r in g_rs.keys():
    a_freqs.append(find_rs_line_array(r, rows))

#a_freqs = ['0' for a in a_freqs if a == 'NaN']
for i,a in enumerate(a_freqs):
    if a == 'NaN':
        a_freqs[i] = 0


a = np.array(a_freqs)
b = a.astype(np.float)

rows_30 =filter(lambda x: getattr(x, 'sample id') =='Mark_2.30', rows)


#!!!!!!WRONG
def get_samples(table, gtype = 'kg'):
    if gtype == 'kg':
        samples = table.columns.keys()
    if gtype == 'hapmap':
        samples = table.columns.keys()
    return samples


l_25 = np.linalg.lstsq(g,b)

a_freqs = []
for r in g_rs.keys():
    a_freqs.append(find_rs_line_array(r, rows_30))

for i,a in enumerate(a_freqs):
    if a == None or a == 'NaN':
        a_freqs[i] = 0

a = np.array(a_freqs)
b = a.astype(np.float)

l_230 = np.linalg.lstsq(g,b)




#kg_samples = get_samples(table = kg_table, gtype = 'kg')


'''
sample_source = {}
for ps in pool_samples:
    if ps.lower() in kg_samples:
        sample_source[ps.lower()] = kg_rsids_sorted
    else:
        print ps
    if ps.lower() not in sample_source.keys() and ps.lower() in hapmap_samples:
        sample_source[ps.lower()] = hapmap_rsids_sorted
'''

'''
def get_genotypes_vcf(sample, rows, index):
    row = rows[index]
    print row
    print sample
    info = getattr(row, sample)
    g = info.split(':')[0]
    g_split = re.split('[/|\\\.|]', g)
    print g_split
    if len(g_split) == 2:
        genotype = sum(map(lambda x: int(x), g_split))
        return genotype
    else:
        return None

def get_all_genotypes(sample_source, index):
    genotypes = []
    for s in sample_source.keys():
        print sample_source[s]
        genotypes.append(get_genotypes_vcf(sample = s, rows = sample_source[s], index = index))
    return genotypes

select_snps = sorted(list(select_rsids))

for i in range(0, len(select_snps)):
    get_all_genotypes(sample_source, index = i)



#s = select([kg_table, hapmap_table], (getattr(kg_table.c, 'rs#') == getattr(hapmap_table.c, 'rs#')), use_labels = True)
#rs = s.execute()
#a_set = set(array_rsids[0:100000])
#a = filter(lambda x: getattr(x, 'kg_lowcov_rs#') in a_set, rs)
#s = select([kg_table, hapmap_table], (getattr(kg_table.c, 'rs#') == getattr(hapmap_table.c, 'rs#')), use_labels = True)


#a = filter(lambda x: x in kg_rsids, list(hapmap_rsids))


#s = session.query(kg_table).filter(getattr(kg_table.c, 'rs#') == inall[1])

''' 