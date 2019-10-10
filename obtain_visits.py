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

def main(tract, patchi, patchj, suffix):
    query = "select DISTINCT(visit), filter from overlaps WHERE tract={} and patch='({}, {})' order by visit".format(tract, patchi, patchj)

    visitab = pd.read_sql_query(query, conn)


    cmd = """nice -n 10 coaddDriver.py /global/cscratch1/sd/bos0109/run2.1i_softln/rerun/calexp-v1 --output $SCRATCH/templates_004 --configfile $HOME/run2_diaproc/coadd_config_example.py  --id tract={} patch={},{} filter={} --selectId visit={} --job templ_03  --cores 4 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-output 
    """
    commands = []
    for filtr, visits in visitab.groupby('filter'):
        visitstr = ''
        for avisit in visits.visit:
            visitstr+=str(avisit)+'^'
            
        commands.append(cmd.format(tract, patchi, patchj, filtr, visitstr[:-1]))

    with open('coaddCommands_{}.sh'.format(suffix), 'w') as cf:
        for acmd in commands:
            cf.write(acmd)
            cf.write('\n \n')
           

if __name__=='__main__':
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
