#array_db.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import commands
import os
import csv 

db = create_engine('sqlite:///../cancergenomes/GENOTYPES.db', echo = True)

Session = sessionmaker(db)
session = Session()
metadata = MetaData(db)

kg_table_name = "1kg_lowcov_1.0"
hapmap_table_name = "hapmap_raw_1.0"

kg_table = Table(kg_table_name, metadata, autoload = True)
hapmap_table = Table(hapmap_table_name, metadata, autoload = True)

afile = open('mark/arraypool/25M1.1', 'r')
lines = afile.readlines()
afile.close()
asnps = set([])
for l in lines:
    asnps.add(l.split('\t')[9] + ':' + l.split('\t')[10])

json.dump(list(asnps), open('array_snps', 'w'))

