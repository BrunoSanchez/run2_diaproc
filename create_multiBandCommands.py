#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  create_multiBandCommands.py
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

repo = "$SCRATCH/templates_rect"
config_path = "./config/multibandDriver_config.py"

nice = "nice -n 10 "
cmd_tmpl = "multiBandDriver.py {} --rerun {} "
cmd_tmpl+= "--id tract={} patch={},{} filter={} "
cmd_tmpl+= "-C {} --config  measureCoaddSources.doPropagateFlags=False "
cmd_tmpl+= "--job {}  --cores {} --time {}  --batch-type={} "
cmd_opt_slrm = " --clobber-config --batch-verbose --batch-stats "
cmd_opt_slrm+= "--mpiexec='-bind-to socket' "
slrm_hasw = "--batch-options='-C haswell -q shared' "
slrm_knl = "--batch-options='-C knl -q regular' "


def main(tract, patch, filters=('g', 'r', 'i', 'z'),
         repo=repo, rerun='multiband', config_path=config_path, 
         cores=None, time=4000, batch='smp', queue_knl=False, 
         outfile=None):
    patchx, patchy = patch 
    strpatch = "'"+str((int(patchx), int(patchy)))+"'"

    fltrstr = ''
    for afiltr in filters:
        fltrstr+=afiltr+'^'
    fltrstr = fltrstr[:-1]
    
    job_name = 'mband_t{}_p{}{}'.format(tract, patchx, patchy)
    if cores is None: cores = len(filters)
    cmd = cmd_tmpl.format(repo, rerun, tract, patchx, patchy, 
                            fltrstr, config_path, job_name, 
                            cores, time, batch)
    
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
            print("Wrote multiband command to {}".format(outfile))
        else:
            try:
                outfile.write(cmd)
                outfile.write('\n\n')    
                print("Wrote multiband command")
            except:
                return cmd
    else:
        return cmd 

    return 


if __name__=='__main__':
    import argparse
    DESC = "Creates commands for multiband merge the coadd images dia_pipe context"
    EPIL = "This will produce a file output with the commands separated by one line"
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)
    parser.add_argument('-t', '--tract', metavar='t', type=int, 
                        help='Tract number')
    parser.add_argument('-p', '--patch', metavar='p', type=str, 
                        help='Patch code number')
    parser.add_argument('-f', '--filters', metavar='f', type=str, 
                        dest='filters', default=['g'])
    parser.add_argument('-c', '--cores', metavar='c', default=None, 
                        help='Number of cores')
    parser.add_argument('-tm', '--time', metavar='tm', default=3000, 
                        help='time limit')
    parser.add_argument('-b', '--batch-type', metavar='b', type=str, 
                        default='smp', help='slurm or smp batch processing')
    parser.add_argument('-q', '--queue-knl', metavar='q', type=str, 
                        default=False, help='true=knl-regular, false=haswell-shared')
    parser.add_argument('-r', '--repo', metavar='r', type=str, default=repo, 
                        help='Repository where templates are located')
    parser.add_argument('-R', '--rerun', metavar='R', type=str, default='multiband', 
                        help='name of rerun for multiband') 
    parser.add_argument('-C', '--conf', metavar='C', type=str, default=config_path, 
                        help='Path of configuration file')
    parser.add_argument('-o', '--outfile', metavar='O', default=None, 
                        help='Outfile path, command output file')
    
    args = parser.parse_args()

    main(args.tract, args.patch, filters=args.filters,
         repo=repo, rerun=args.rerun, config_path=args.conf, 
         cores=args.cores, time=args.time, batch=args.batch_type, 
         queue_knl=args.queue_knl, outfile=args.outfile)