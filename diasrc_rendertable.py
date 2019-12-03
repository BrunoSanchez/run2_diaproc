#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  diasrc_rendertable.py
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
# NEEDS TO RUN USING desc-dia 

import os
import sqlite3

import numpy as np

import lsst.afw.cameraGeom
import lsst.afw.geom as afwGeom
import lsst.geom as geom

import matplotlib.pyplot as plt
import pandas as pd

from astropy import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import vstack

from lsst.afw.geom import makeSkyWcs
from lsst.daf.persistence import Butler
from lsst.obs.lsst.imsim import ImsimMapper

calexprepo = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1' 
b = Butler(calexprepo)
skymap = b.get('deepCoadd_skyMap')

template_repo = '/global/cscratch1/sd/bos0109/templates_rect'
diarepo = template_repo + '/rerun/diff_rect'
assocrepo = diarepo + '/rerun/assoc_sha'
forcerepo = assocrepo + '/rerun/forcedPhot' 
tmprepo = template_repo + '/rerun/multiband'

diabutler = Butler(forcerepo)

#truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58.0_56.0_-31.0_-32.0.csv')
#sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58.0_56.0_-31.0_-32.0.csv')

diaSrc_store = pd.HDFStore('/global/cscratch1/sd/bos0109/diaSrc_fulltables_v2.h5')
diaSrc_store.open()
metacols = ['id', 'visit', 'filter', 'raftName', 'detectorName', 'detector']

## ========================================================================== ##
## =================== Matching iterating over t+p ========================== ##
## ========================================================================== ##

## First a rectangle of (alpha, delta)
ramax = 58
ramin = 56
decmax = -31
decmin = -32

## from here we can pick the list of tract+patches
# creating a rectangle of 2 sq. degree for tract/patch search
radec_NE = afwGeom.SpherePoint(ramax, decmax, afwGeom.degrees)
radec_SE = afwGeom.SpherePoint(ramax, decmin, afwGeom.degrees)
radec_SW = afwGeom.SpherePoint(ramin, decmin, afwGeom.degrees)
radec_NW = afwGeom.SpherePoint(ramin, decmax, afwGeom.degrees)
rect = [radec_NE, radec_NW, radec_SW, radec_SE]

tpatches = skymap.findTractPatchList(rect)

# iterating over tract and patches to build diaSrc catalogues
path = os.path.join(diarepo, 'deepDiff')
diffpath = 'v{}-f{}/{}/diaSrc_{}-{}-{}-{}-det{}.fits'
for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}
        
        store_key = str(tract.getId())+'_'+str(patch_i)+str(patch_j)
        if store_key in diaSrc_store: continue
        
        print('starting with ', tpId)
        metadata = diabutler.queryMetadata('deepDiff_diaSrc',metacols,dataId=tpId)
        metadata = pd.DataFrame(metadata, columns=metacols)
        #metadata = metadata[metadata['filter']!='u']
        metadata = metadata[metadata['filter']!='y']
        
        cats = []
        for idx, idr, vn, fna, raf, detN, det in metadata.itertuples():
            #if fna=='y' or fna=='u': continue
            pp = diffpath.format(str(vn).zfill(8), fna, raf, str(vn).zfill(8), 
                                 fna, raf, detN, str(det).zfill(3))
            dpath = os.path.join(path, pp)
            if os.path.exists(dpath):
                catalog = diabutler.get('deepDiff_diaSrc', visit=vn, 
                                        detector=det).asAstropy()
                if len(catalog) is not 0:
                    catalog['visit_n'] = vn
                    catalog['filter'] = fna
                    catalog['raft'] = raf
                    catalog['sensor'] = detN
                    catalog['detector'] = det
                    cats.append(catalog)
        mastercat = vstack(cats)
        mastercat = mastercat.to_pandas()
        diaSrc_store[store_key] = mastercat
        diaSrc_store.flush()
diaSrc_store.close()
