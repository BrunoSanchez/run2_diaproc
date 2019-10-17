
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

from glob import glob
import pandas as pd

repo_dir = '/global/cscratch1/sd/bos0109/test_imdiff_run2/deepDiff/'

templ = 'v*'


cmd = "time nice -n 10 forcedPhotCcdDiaDriver.py /global/cscratch1/sd/bos0109/test_imdiff_run2/rerun/multiband --rerun forcedPhot  --id visit={} --cores 4"

visits_str = ''
for adir in glob(repo_dir+templ):
    slimdir = adir.split('/')[-1]
    visitn = slimdir[1:-3]
    print(slimdir, visitn)
    visits_str += visitn+'^'

commands = [cmd.format(visits_str[:-1])]

with open('forcedPhotDiaCommands.sh', 'w') as cf:
    for acmd in commands:
        cf.write(acmd)
        cf.write('\n \n')
        
