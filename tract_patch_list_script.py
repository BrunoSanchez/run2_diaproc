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

# check the tracts+patch, print their names
for atract in tpatches:
    print(atract[0])  # prints the number of the tract
    for thepatches in atract[1:]:  # next things on list are the patches
        for apatch in thepatches:
            print(apatch)
            # Now we have the list of tract+patch 
            # let's find the visit list
            patchx, patchy = apatch.getIndex()
            strpatch = "'"+str((int(patchx), int(patchy)))+"'"
            query = query_tmpl.format(atract[0].getId(), strpatch)
            visitab = pd.read_sql_query(query, conn)
            
            ccoadd.main(atract[0].getId(), 
                       apatch.getIndex(), calexp_repo=repo,
                       output_repo="$SCRATCH/templates_005", 
                       database=database, cores=4, batch='slurm')
            

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

