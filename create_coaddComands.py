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
query_tmpl = "select DISTINCT(visit), filter from overlaps WHERE tract={} and patch={} order by visit"
calexp_repo = "/global/cscratch1/sd/bos0109/run2.1i_softln/rerun/calexp-v1"
output_repo = "$SCRATCH/templates_004"
config_path = "$HOME/run2_diaproc/coadd_config_example.py"

cmd_tmpl = "nice -n 10 coaddDriver.py {} --output {} --configfile {} "
cmd_tmpl+= "--id tract={} patch={},{} filter={} --selectId visit={} "
cmd_tmpl+= "--job {}  --cores {} --time 600  --batch-type={} "
cmd_tmpl+= "#  --batch-verbose  --batch-stats "
cmd_tmpl+= "--batch-options='-C knl -q regular' --mpiexec='-bind-to socket' "
cmd_tmpl+= "  --clobber-output --dry-run"

def main(tract, patch, calexp_repo=calexp_repo,
         output_repo=output_repo, config_path=config_path, 
         database=database, cores=4, batch='smp'):
    conn = sqlite3.connect(database)
    c = conn.cursor()
    patchx, _, patchy = patch 
    strpatch = "'"+str((int(patchx), int(patchy)))+"'"
    query = query_tmpl.format(tract, strpatch)
    visitab = pd.read_sql_query(query, conn)

    commands = []
    for filtr, visits in visitab.groupby('filter'):
        visitstr = ''
        for avisit in visits.visit:
            visitstr+=str(avisit)+'^'
        job_name = 'coadd_t{}_p{}{}'.format(tract, patchx, patchy)
        cmd = cmd_tmpl.format(calexp_repo, output_repo, config_path,
            tract, patchx, patchy, filtr, visitstr[:-1], job_name, 
            cores, batch)
        commands.append(cmd)
    outfile = 'driver_commands/coaddCommands_t{}_p{}{}.sh'.format(
        tract, patchx, patchy)
    with open(outfile, 'w') as cf:
        for acmd in commands:
            cf.write(acmd)
            cf.write('\n \n')
    print("Wrote {} commands to {}".format(len(commands), outfile))
    return      


if __name__=='__main__':
    import argparse
    DESC = "Creates commands for coadd image driver in the dia_pipe context"
    EPIL = "This will produce a file output with the commands separated by one line"
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)
    parser.add_argument('-t', '--tract', metavar='t', type=int, 
                        help='Tract number')
    parser.add_argument('-p', '--patch', metavar='p', type=str, 
                        help='Patch code number')
    parser.add_argument('-db', '--database', metavar='db', type=str, 
                        dest='database', default=database, 
                        help='Database of tract+patchs to find out visits')
    parser.add_argument('-c', '--cores', metavar='c', type=int, default=4, 
                        help='Number of cores')
    parser.add_argument('-b', '--batch-type', metavar='b', type=str, 
                        default='smp', help='slurm or smp batch processing')
    parser.add_argument('-x', '--calexp', metavar='Cx', type=str, default=calexp_repo, 
                        help='Repository where calexps are located')
    parser.add_argument('-r', '--tmpl', metavar='R', type=str, default=output_repo, 
                        help='Repository where templates are going to be located')
    parser.add_argument('-C', '--conf', metavar='C', type=str, default=config_path, 
                        help='Path of configuration file')
    args = parser.parse_args()

    main(args.tract, args.patch, database=args.database, 
         batch=args.batch_type, cores=args.cores, 
         calexp_repo=args.calexp, output_repo=args.tmpl, 
         config_path=args.conf)
