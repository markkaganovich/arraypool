#array_db.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table, select
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

def get_rsids(name, table):
    if name in os.listdir('./'):
        print "loading from file"
        rsids = json.load(open(name, 'r'))
    else:
        rsids = set([])
        s = table.select()
        rs = s.execute()
        for row in rs:
            rsids.add(getattr(row, 'rs#'))
        json.dump(list(rsids), open(name, 'w'))
    return rsids

kg_table = Table(kg_table_name, metadata, autoload = True)
hapmap_table = Table(hapmap_table_name, metadata, autoload = True)
array_table = Table(array_table_name, metadata_array, autoload = True)

print "post-table"
#inboth = kg_rsids.intersection(hapmap_rsids)

s = select([kg_table, hapmap_table, array_table], (getattr(kg_table.c, 'rs#') == getattr(hapmap_table.c, 'rs#')) and (getattr(kg_table.c, 'rs#') == getattr(array_table.c, 'snp name')) )
rs = s.execute()

#a = filter(lambda x: x in kg_rsids, list(hapmap_rsids))


