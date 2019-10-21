
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
query_tmpl = "select DISTINCT(visit), filter from overlaps WHERE tract={} and patch={} order by visit"
cmd = "time nice -n 10 imageDifferenceDriver.py  /global/cscratch1/sd/bos0109/templates_003/rerun/multiband \
       --output /global/cscratch1/sd/bos0109/test_imdiff_run2  --id visit={} \
       -C imageDifferenceDriver_config.py --batch-type={} --mpiexec='-bind-to socket'\
       --cores {}  --job imdiff_v{}_f{} --time 500 --batch-options='-C knl -q regular'"


def main(tract, patch, filters='griz', 
         outfile='driver_commands/diaCommands.sh', 
         database=database, batch='smp', cores=4):
    conn = sqlite3.connect(database)
    c = conn.cursor()
    visitab = pd.read_sql_query(query.format(tract, str(patch)), conn)
    commands = []
    for filtr, visits in visitab.groupby('filter'):
        if filtr in list(filters):
            print(filtr, visits.visit)
            for avisit in visits.visit:
                commands.append(cmd.format(avisit, batch, cores, avisit, filtr))

    with open(outfile, 'w') as cf:
        for acmd in commands:
            cf.write(acmd)
            cf.write('\n \n')
    print("Wrote {} commands to {}".format(len(commands), outfile))
    return



if __name__=='__main__':
    import sys
    import argparse
    DESC = """Creates commands for difference image driver in the dia_pipe context."""
    EPIL = """This will produce a file output with the commands separated by one line."""
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)
    parser.add_argument('-t', '--tract', metavar='t', type=int, 
                        help='Tract number')
    parser.add_argument('-p', '--patch', metavar='p', type=str, 
                        help='Patch code number')
    parser.add_argument('-f', '--filter', metavar='f', type=str, help='Filter name')
    parser.add_argument('-o', '--outfile', metavar='o', type=str, help='Output filename with the commands')
    parser.add_argument('-db', '--database', metavar='db', type=str, dest='database', 
                        default=database, help='Database of tract+patchs to find out visits')
    parser.add_argument('-c', '--cores', metavar='c', type=int, default=4, help='Number of cores')
    parser.add_argument('-b', '--batch-type', metavar='b', type=str, default='smp', help='slurm or smp batch processing')

    args = parser.parse_args()
    print(args)
    main(args.tract, args.patch, filters=args.filter, 
         outfile=args.outfile, database=args.database, 
         batch=args.batch_type, cores=args.cores)
