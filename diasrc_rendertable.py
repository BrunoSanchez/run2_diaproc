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

#import lsst.afw.cameraGeom
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
tmprepo = template_repo + '/rerun/multiband'

diarepo = template_repo + '/rerun/diff_rect'
assocrepo = diarepo + '/rerun/assoc_thirrun'
forcerepo = assocrepo + '/rerun/forcedPhot' 

diabutler = Butler(forcerepo)

#truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58.0_56.0_-31.0_-32.0.csv')
#sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58.0_56.0_-31.0_-32.0.csv')

diaSrc_store = pd.HDFStore('/global/homes/b/bos0109/run2_diaproc/results/diaSrc_secrun_fulltables_v4.h5')
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

# region ----------------------------------------------------------------------------------------------
# iterating over tract and patches to build diaForced catalogues
metas = []
for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}
      
        #store_key = str(tract.getId())+'_'+str(patch_i)+str(patch_j)
        #if store_key in diaSrc_store: continue
        print('starting with ', tpId)
        metadata = diabutler.queryMetadata('deepDiff_forced_diaSrc', metacols,dataId=tpId)
        metadata = pd.DataFrame(metadata, columns=metacols)
        metadata = metadata[metadata['filter']!='u']
        metadata = metadata[metadata['filter']!='y']
        #metadata['tract'] = tpId['tract']
        #metadata['patch'] = tpId['patch']
        metas.append(metadata)
metadata = pd.concat(metas).drop_duplicates()
cats = []
path = os.path.join(forcerepo, 'deepDiff')
#diffpath = 'v{}-f{}/{}/diaSrc_{}-{}-{}-{}-det{}.fits'
diffpath = 'v{}-f{}/{}/diaForced_{}-{}-{}-{}-det{}.fits'
for idx, idr, vn, fna, raf, detN, det, in metadata.itertuples():
    #if fna=='y' or fna=='u': continue
    pp = diffpath.format(str(vn).zfill(8), fna, raf, str(vn).zfill(8), 
                            fna, raf, detN, str(det).zfill(3))
    dpath = os.path.join(path, pp)
    if os.path.exists(dpath):
        catalog = diabutler.get('deepDiff_forced_diaSrc', visit=vn, 
                #tract=int(t), patch=p, 
                detector=det).asAstropy()
        if len(catalog) != 0:
            catalog['visit_n'] = vn
            catalog['filter'] = fna
            catalog['raft'] = raf
            catalog['sensor'] = detN
            catalog['detector'] = det
            cats.append(catalog)
#import ipdb; ipdb.set_trace()
mastercat = vstack(cats)
mastercat = mastercat.to_pandas()
#diaSrc_store[store_key] = mastercat
diaSrc_store['full_table_forced'] = mastercat
diaSrc_store.flush()
# endregion -------------------------------------------------------------------------------------------
del(mastercat)
del(cats)
# iterating over tract and patches to build diaSrc catalogues
# region ----------------------------------------------------------------------------------------------
metas = []
for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}
        
        #store_key = str(tract.getId())+'_'+str(patch_i)+str(patch_j)
        #if store_key in diaSrc_store: continue
        print('starting with ', tpId)
        metadata = diabutler.queryMetadata('deepDiff_diaSrc', metacols,dataId=tpId)
        metadata = pd.DataFrame(metadata, columns=metacols)
        metadata = metadata[metadata['filter']!='u']
        metadata = metadata[metadata['filter']!='y']
        #metadata['tract'] = tpId['tract']
        #metadata['patch'] = tpId['patch']
        metas.append(metadata)
metadata = pd.concat(metas).drop_duplicates()

cats = []
path = os.path.join(diarepo, 'deepDiff')
diffpath = 'v{}-f{}/{}/diaSrc_{}-{}-{}-{}-det{}.fits'
#diffpath = 'v{}-f{}/{}/diaForced_{}-{}-{}-{}-det{}.fits'
for idx, idr, vn, fna, raf, detN, det in metadata.itertuples():
    #if fna=='y' or fna=='u': continue
    pp = diffpath.format(str(vn).zfill(8), fna, raf, str(vn).zfill(8), 
                         fna, raf, detN, str(det).zfill(3))
    dpath = os.path.join(path, pp)
    if os.path.exists(dpath):
        catalog = diabutler.get('deepDiff_diaSrc', visit=vn, 
                #tract=int(t), patch=p, 
                detector=det).asAstropy()
        if len(catalog) != 0:
            catalog['visit_n'] = vn
            catalog['filter'] = fna
            catalog['raft'] = raf
            catalog['sensor'] = detN
            catalog['detector'] = det
            cats.append(catalog)
#import ipdb; ipdb.set_trace()
mastercat = vstack(cats)
mastercat = mastercat.to_pandas()
#diaSrc_store[store_key] = mastercat
diaSrc_store['full_table'] = mastercat
diaSrc_store.flush()
# endregion -------------------------------------------------------------------------------------------
diaSrc_store.close()