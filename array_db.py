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

print "here"

kg_table = Table(kg_table_name, metadata, autoload = True)
hapmap_table = Table(hapmap_table_name, metadata, autoload = True)
#array_table = Table(array_table_name, metadata_array, autoload = True)

print "post-table"

if 'hapmap_rsids' in os.listdir('./'):
    print "loading from file"
    hapmap_rsids = set(json.load(open('hapmap_rsids')))
else:
    hapmap_rsids = set([])
    hm = hapmap_table.select()
    rs = hm.execute()
    for row in rs:
        hapmap_rsids.add(getattr(row, 'rs#'))
    json.dump(list(hapmap_rsids), open('hapmap_rsids', 'w'))

if 'kg_rsids' in os.listdir('./'):
    print "loading from file"
    kg_rsids = set(json.load(open('kg_rsids')))
else:
    kg_rsids = set([])
    kg = kg_table.select()
    rs = kg.execute()
    for row in rs:
        kg_rsids.add(getattr(row, 'rs#'))
    json.dump(list(kg_rsids), open('kg_rsids', 'w'))

inboth = kg_rsids.intersect(hg_rsids)

#a = filter(lambda x: x in kg_rsids, list(hapmap_rsids))


