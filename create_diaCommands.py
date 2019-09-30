
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  create_diaCommands.py
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

cmd = "time nice -n 10 imageDifferenceDriver.py  /global/cscratch1/sd/bos0109/templates_003/rerun/multiband --output /global/cscratch1/sd/bos0109/test_imdiff_run2  --id visit={} -C imageDifferenceDriver_config.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 4  --job imdiff_v{}_f{} --time 500 --batch-options='-C knl -q regular'"


commands = []
for filtr, visits in visitab.groupby('filter'):
    if filtr not in ['u','y']:
        print(filtr, visits.visit)
        for avisit in visits.visit:
            commands.append(cmd.format(avisit, avisit, filtr))

with open('diaCommands.sh', 'w') as cf:
    for acmd in commands:
        cf.write(acmd)
        cf.write('\n \n')
        
