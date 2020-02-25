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
import pandas as pd

rerun = 'assoc'
diff_repo = '$SCRATCH/templates_rect/rerun/diff_rect'

nice = "time nice -n 10 "
cmd_tmpl = "associationDriver.py  {} --rerun {} "
cmd_tmpl +="--id tract={} filter={} --selectId visit={} --batch-type={} " 
cmd_tmpl +="--clobber-config --cores {} --job assoc{} --time {} "
cmd_opt_slrm = " --batch-verbose  --batch-stats "
cmd_opt_slrm+= "--mpiexec='-bind-to socket' "
slrm_hasw = "--batch-options='-C haswell -q shared' "
slrm_knl = "--batch-options='-C knl -q regular' "


def main(tract, filters='griz', diff_repo=diff_repo, visitab=None,
        outfile='driver_commands/assocCommands.sh', 
        batch='smp', cores=None, rerun=rerun, queue_knl=False, time=100):

    visits_str = ''
    #file_prefix_templ = 'v*'
    if visitab is None:
        for adir in glob(diff_repo+'/deepDiff/v*'):
            slimdir = adir.split('/')[-1]
            visitn = slimdir[1:-3]
            print(slimdir, visitn)
            visits_str += visitn+'^'
    elif isinstance(visitab, pd.DataFrame):
        #import ipdb; ipdb.set_trace()
        for visitn in visitab.visit.values:
            visits_str += str(visitn)+'^'

    fltrstr = ''
    for afiltr in list(filters):
        fltrstr+=afiltr+'^'
    fltrstr = fltrstr[:-1]
    if cores is None: cores = 8

    cmd = cmd_tmpl.format(diff_repo, rerun, tract, fltrstr, 
                          visits_str[:-1], batch, cores, tract, time)

    if batch=='slurm':
        cmd+=cmd_opt_slrm
        if queue_knl:
            cmd+=slrm_knl
        else:
            cmd+=slrm_hasw
    else: 
        cmd = nice + cmd
    
    if outfile is not None:
        if isinstance(outfile, str):
            with open(outfile, 'a+') as cf:
                cf.write(cmd)
                cf.write('\n\n')
            print("Wrote association command to {}".format(outfile))
        else:
            try:
                outfile.write(cmd)
                outfile.write('\n\n')    
                print("Wrote association command")
            except:
                return cmd

    return


if __name__=='__main__':
    import argparse
    DESC = "Creates commands for association driver in the dia_pipe context"
    EPIL = "This will produce a file output with the commands separated by one line"
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)
    parser.add_argument('-t', '--tract', metavar='t', type=int, 
                        help='Tract number')
    parser.add_argument('-f', '--filter', metavar='f', type=str, 
                        help='Filter names', default='griz')
    parser.add_argument('-o', '--outfile', metavar='o', type=str, 
                        help='Output filename with the commands', 
                        default='driver_commands/diaCommands.sh')
    parser.add_argument('-c', '--cores', metavar='c', type=int, default=4, 
                        help='Number of cores')
    parser.add_argument('-b', '--batch-type', metavar='b', type=str, 
                        default='smp', help='slurm or smp batch processing')
    parser.add_argument('-d', '--diff', metavar='d', type=str, default=diff_repo, 
                        help='Repository where templates are located')
    parser.add_argument('-r', '--rerun', metavar='re', type=str, default=rerun, 
                        help='Repository where assoc cats are going to be located')
    parser.add_argument('-tm', '--time', metavar='tm', type=int, 
                        default=120, help='time of execution per visit')
    parser.add_argument('-q', '--queue', metavar='q', default=False, 
                        help='True: knl -q regular; False: haswell -q shared')
                        
    args = parser.parse_args()

    main(args.tract, filters=args.filter, outfile=args.outfile, 
         batch=args.batch_type, cores=args.cores, 
         diff_repo=args.diff_repo, rerun=args.rerun, 
         queue_knl=args.queue, time=args.time)
