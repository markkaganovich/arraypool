#array_db.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table, select
import sys
import json
import commands
import os
import csv 

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

kg_table_name = "kg"
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

ceu = KG('ceu', kg_files[0], 0)
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
inall = list(inboth.intersection(set(inboth)))
select_rsids = set(inall[1:1000])



def get_rows(attr = 'rs#', selected_attr = [], table = None):
    rows = []
    s = table.select()
    rs = s.execute()
    for r in rs:
        if getattr(r, attr) in selected_attr:
            rows.append(r)
    return rows

#kg_rsid_rows = get_rows(attr = 'rs#', selected_attr = select_rsids, table = kg_table)
hapmap_rsid_rows = get_rows(attr = 'rs#', selected_attr = select_rsids, table = hapmap_table)

#kg_rsids_sorted = sorted(kg_rsid_rows, key=lambda r: getattr(r, 'rs#'))
hapmap_rsids_rows_sorted = sorted(hapmap_rsid_rows, key=lambda r: getattr(r, 'rs#'))

for i in range(1,len(hapmap_rsid_rows)-1): 
    if getattr(hapmap_rsids_sorted[i-1], 'rs#') == getattr(hapmap_rsids_sorted[i], 'rs#') \
        or (getattr(hapmap_rsids_sorted[i+1], 'rs#') == getattr(hapmap_rsids_sorted[i], 'rs#')):
        hapmap_rsid_rows.remove(hapmap_rsids_sorted[i])
'''
given list of cell line IDs, find if they are coming from hapmap or 1kg, 
get their GENOTYPES for all the rsids 
'''
pool1_samples = ["NA18516", "NA18517", "NA18579", "NA18592", "NA18561", "NA07357", "NA06994", "NA18526", "NA12004", "NA19141", "NA19143", "NA19147", "NA19152", "NA19153", "NA19159", "NA19171", "NA19172", "NA19190", "NA19207", "NA19209", "NA19210", "NA19225", "NA18856", "NA18858", "NA18562", "NA18563", "NA18853", "NA18861"]

pool_samples = pool1_samples
hapmap_samples = get_samples(table = hapmap_table, gtype = 'hapmap')


def find_rs_line_kg(k_object, rs, s):
    for l in k_object.lines:
        if rs == l.split('\t')[k_object.header.index('ID')]:                
            info = l[k_object.header.index(s)]
            g = info.split(':')[0]
            g_split = re.split('[/|\\\.|]', g)
            if len(g_split) == 2:
                genotype = sum(map(lambda x: int(x), g_split))
                return genotype

    '''            
    if source == 'hapmap':
        for h in hapmap_rsids_rows_sorted:
            if getattr(h, 'rs#') == rs:
    '''



rs = select_snps[0]
rs_lines = []
for p in pool_samples:
    for k in kgs:
        if p in k.header:
            source = 'kg'
            find_rs_line_kg(k, rs, p)         
        #if p in hapmap_samples:
        #    source = 'hapmap'
        #    find_rs_line()
        #print source

#!!!!!!WRONG
def get_samples(table, gtype = 'kg'):
    if gtype == 'kg':
        samples = table.columns.keys()
    if gtype == 'hapmap':
        samples = table.columns.keys()
    return samples

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