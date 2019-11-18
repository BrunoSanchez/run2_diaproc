
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# create_forcedPhotCommands.py
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

from glob import glob
import pandas as pd

dia_repo_parent = '/global/cscratch1/sd/bos0109/test_imdiff_run2/deepDiff/'

dia_repo = "/global/cscratch1/sd/bos0109/test_imdiff_run2/rerun/multiband"

cmd_tmpl = "time nice -n 10 forcedPhotCcdDiaDriver.py {} --rerun forcedPhot "
cmd_tmpl +=" --id visit={} --cores {} --batch-type={} "
cmd_tmpl +="--batch-options='-C knl -q regular'"

def main(dia_repo=dia_repo, dia_parent=dia_repo_parent, 
         outfile='forcedPhotDiaCommands.sh', 
         cores=4, batch_type='smp'):
    visits_str = ''
    #file_prefix_templ = 'v*'
    for adir in glob(dia_repo_parent+'v*'):
        slimdir = adir.split('/')[-1]
        visitn = slimdir[1:-3]
        print(slimdir, visitn)
        visits_str += visitn+'^'

    commands = [cmd_tmpl.format(dia_repo, visits_str[:-1], cores, batch_type)]

    with open(outfile, 'w') as cf:
        for acmd in commands:
            cf.write(acmd)
            cf.write('\n \n')
    
    print("Wrote {} commands to {}".format(len(commands), outfile))
    return
    

if __name__=='__main__':

    import argparse
    DESC = "Creates commands for forced photometry driver in the dia_pipe context"
    EPIL = "This will produce a file output with the commands separated by one line"
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)
    parser.add_argument('-o', '--outfile', metavar='o', type=str, 
                        help='Output filename with the commands', 
                        default='driver_commands/forcedPhotDiaCommands.sh')
    parser.add_argument('-c', '--cores', metavar='c', type=int, default=4, 
                        help='Number of cores')
    parser.add_argument('-b', '--batch-type', metavar='b', type=str, 
                        default='smp', help='slurm or smp batch processing')
    parser.add_argument('-d', '--diff', metavar='D',
                        type=str, default=dia_repo_parent, 
                        help='Repository where dia first files are located')
    parser.add_argument('-r', '--rerun', metavar='R', type=str, default=dia_repo, 
                        help='Repository where dia last rerun is located')
    args = parser.parse_args()

    main(dia_repo=args.rerun, dia_parent=args.diff, 
         outfile=args.outfile, cores=args.cores, batch_type=args.batch_type)
