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

afile = open('mark/arraypool/25M1.1', 'r')
lines = afile.readlines()
afile.close()
asnps = set([])
for l in lines:
    asnps.add(l.split('\t')[9] + ':' + l.split('\t')[10])

json.dump(asnps, open('array_snps', 'w'))

