#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  obtain_visits.py
#  
#  Copyright 2019 bruno <bruno.sanchez@duke.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sqlite3
import pandas as pd

database = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1/tracts_mapping.sqlite3'

conn = sqlite3.connect(database)
c = conn.cursor()

query = "select DISTINCT(visit), filter from overlaps WHERE tract=4639 and patch='(0, 0)' order by visit"

visitab = pd.read_sql_query(query, conn)


cmd = """nice -n 10 coaddDriver.py /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1 --output $SCRATCH/templates_003 --configfile $HOME/run2_diaproc/coadd_config_example.py  --id tract=4639 patch=0,0 filter={} --selectId visit={} --job templ_02  --cores 32 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-output 
"""
commands = []
for filtr, visits in visitab.groupby('filter'):
    visitstr = ''
    for avisit in visits.visit:
        visitstr+=str(avisit)+'^'
        
    commands.append(cmd.format(filtr, visitstr[:-1]))

with open('coaddCommands.sh', 'w') as cf:
    for acmd in commands:
        cf.write(acmd)
        cf.write('\n \n')
        
