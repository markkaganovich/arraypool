#array_db.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import commands
import os
import csv 

db = create_engine('sqlite:///../cancergenomes/GENOTYPES.db', echo = False)
db_array = create_engine('sqlite:///../arraydata/arrays.db')


Session = sessionmaker(db)
session = Session()
metadata = MetaData(db)

Session_array = sessionmaker(db_array)
session_array = Session_array()
metadata_array = MetaData(db_array)

kg_table_name = "kg_lowcov"
hapmap_table_name = "hapmap"
array_table_name = "Array1S_finalreport"

def copy_my_table(name, table):
    args = []
    for c in table.columns:
        args.append(c.copy())
    for c in table.constraints:
        args.append(c.copy())
    return Table(name, table.metadata, *args)



kg_table = Table(kg_table_name, metadata, autoload = True)
hapmap_table = Table(hapmap_table_name, metadata, autoload = True)
array_table = Table(array_table_name, metadata_array, autoload = True)


hapmap_rsids = set([])
for row in session.query(hapmap_table):
    hapmap_rsids.add(getattr(row, 'rs#'))


'''
afile = open('mark/arraypool/25M1.1', 'r')
lines = afile.readlines()
afile.close()
asnps = set([])
for l in lines:
    asnps.add(l.split('\t')[9] + ':' + l.split('\t')[10])

json.dump(list(asnps), open('array_snps', 'w'))
'''

json.dump(list(hapmap_rsids), open('hapmap_rsids', 'w'))

