#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  tract_patch_list_script.py
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

import os
import glob
import warnings
import sqlite3
import re

import numpy as np
import pandas as pd

from astropy.table import Table

import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisplay
import lsst.afw.cameraGeom as cameraGeom

from lsst.daf.persistence import Butler

from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import create_coaddComands as ccoadd 
import create_multiBandCommands as multib
import create_diaCommands as cdia
import create_assocCommands  as cassoc
import create_forcedPhotCommands as cfPhot 

repo = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1' 
b = Butler(repo)

 
skymap = b.get('deepCoadd_skyMap')

# creating a rectangle of 2 sq. degree for tract/patch search
radec_NE = afwGeom.SpherePoint(58, -31, afwGeom.degrees)
radec_SE = afwGeom.SpherePoint(58, -32, afwGeom.degrees)
radec_SW = afwGeom.SpherePoint(56, -32, afwGeom.degrees)
radec_NW = afwGeom.SpherePoint(56, -31, afwGeom.degrees)
rect = [radec_NE, radec_NW, radec_SW, radec_SE]

tpatches = skymap.findTractPatchList(rect)

# to figure the visits on this patches we need the DB  
# this will help us dtermine seeing and time cuts for coadds
database = repo+'/tracts_mapping.sqlite3'
query_tmpl = "select DISTINCT(visit), filter from overlaps WHERE tract={} and patch={} order by visit"
conn = sqlite3.connect(database)
c = conn.cursor()

#coaddOutfile = './driver_commands/coadd_multiband_coadd.sh'
multibOutfile = './driver_commands/coadd_multiband_coadd.sh'
diaOutfile = 'driver_commands/diaCommands_fullset.sh'
assocOutfile = 'driver_commands/assocCommands_fullset.sh'
forcedOutfile = 'driver_commands/forcedPhotCommands.sh'

for acmdfile in [multibOutfile, diaOutfile, assocOutfile]:
    if os.path.exists(acmdfile):
        os.remove(acmdfile)

full_visits = []
# check the tracts+patch, print their names
for atract in tpatches:
    print(atract[0])  # prints the number of the tract
    thepatches = atract[1]  # next things on list are the patches
    for apatch in thepatches:
        print(apatch)
        # Now we have the list of tract+patch 
        # let's find the visit list
        patchx, patchy = apatch.getIndex()
        strpatch = "'"+str((int(patchx), int(patchy)))+"'"
        query = query_tmpl.format(atract[0].getId(), strpatch)
        visitab = pd.read_sql_query(query, conn)
        full_visits.append(visitab)

        ccoadd.main(atract[0].getId(), apatch.getIndex(), calexp_repo=repo,
                    output_repo="$SCRATCH/templates_rect", 
                    database=database, cores=4, batch='slurm', 
                    queue_knl=True)
        
        multib.main(atract[0].getId(), apatch.getIndex(), 
                    filters='ugriz',
                    repo="$SCRATCH/templates_rect", 
                    rerun='multiband', batch='slurm', queue_knl=True,
                    outfile=multibOutfile)
    cassoc.main(atract[0].getId(), filters='ugriz', outfile=assocOutfile,
                diff_repo="$SCRATCH/templates_rect/rerun/diff_rect",
                batch='smp', cores=len(thepatches), rerun='assoc_sha', time=500)

full_visitab = pd.concat(full_visits).drop_duplicates('visit').reset_index(drop=True)
full_visitab.to_csv('./catalogs+tables/full_visits_from_tractmapping_db.csv')
          
cdia.main(filters='ugriz', visit=full_visitab,
          outfile=diaOutfile, batch='slurm', cores=4, queue_knl=True,
          tmpl_repo="$SCRATCH/templates_rect", rerun="diff_rect", 
          config_path="./config/imageDifferenceDriver_config.py")

cfPhot.main(dia_repo="$SCRATCH/templates_rect/rerun/diff_rect/rerun/assoc_sha",
            dia_parent="$SCRATCH/templates_rect/rerun/diff_rect", time=50,
            outfile=forcedOutfile, cores=8, batch_type='slurm', queue_knl=True)
            
            
            

## I can get the tract and patch for every single patch in the sky area
## now I should feed this to my script to launch several jobs
## to do this I must import my own module


# Mapper = b._getDefaultMapper()
# mapper = Mapper(repo)
# all_diadataset_types = mapper.getDatasetTypes() 
 
# remove = ['_config', '_filename', '_md', '_sub', '_len', '_schema', '_metadata'] 
 
# shortlist = [] 
# for dataset_type in all_diadataset_types: 
#     keep = True 
#     for word in remove: 
#         if word in dataset_type: 
#             keep = False 
#     if keep: 
#         shortlist.append(dataset_type) 
# for adtype in shortlist: 
#     try: 
#         ttypes = b.getKeys(adtype)
#         if len(ttypes) is not 0:
#             print('++++++++++++\n', adtype) 
#             print(ttypes) 
#     except KeyError: 
#         continue  

