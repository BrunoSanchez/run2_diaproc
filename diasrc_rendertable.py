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
os.environ['SCRATCH']='/global/cscratch1/sd/bos0109'
SCRATCH = '/global/cscratch1/sd/bos0109'

import sqlite3

import numpy as np

#import lsst.afw.cameraGeom
#import lsst.afw.geom as afwGeom   ## deprecated
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
store = f'{SCRATCH}/results/diaSrc_thirrun_fulltables_v1.h5'
diaSrc_store = pd.HDFStore(store)
diaSrc_store.open()
metacols = ['id', 'visit', 'filter', 'raftName', 'detectorName', 'detector']

# =============================================================================
# optimizing df
# =============================================================================
# We're going to be calculating memory usage a lot,
# so we'll create a function to save us some time!
def mem_usage(pandas_obj):
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb)


def optimize_df(df):
    df_int = df.select_dtypes(include=['int'])
    #converted_int = df_int.apply(pd.to_numeric, downcast='unsigned')

    #print(mem_usage(df_int))
    #print(mem_usage(converted_int))

    compare_ints = pd.concat([df_int.dtypes,converted_int.dtypes],axis=1)
    compare_ints.columns = ['before','after']
    compare_ints.apply(pd.Series.value_counts)

    df_float = df.select_dtypes(include=['float'])
    converted_float = df_float.apply(pd.to_numeric,downcast='float')

    #print(mem_usage(df_float))
    #print(mem_usage(converted_float))

    compare_floats = pd.concat([df_float.dtypes,converted_float.dtypes],axis=1)
    compare_floats.columns = ['before','after']
    compare_floats.apply(pd.Series.value_counts)

    optimized_df = df.copy()

    optimized_df[converted_int.columns] = converted_int
    optimized_df[converted_float.columns] = converted_float

    mem_df = mem_usage(df)
    mem_op_df = mem_usage(optimized_df)
    print(mem_df)
    print(mem_op_df)
    if mem_df<=mem_op_df:
        print('Memory increased, returning original')
        return df

    return optimized_df

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
radec_NE = geom.SpherePoint(ramax, decmax, geom.degrees)
radec_SE = geom.SpherePoint(ramax, decmin, geom.degrees)
radec_SW = geom.SpherePoint(ramin, decmin, geom.degrees)
radec_NW = geom.SpherePoint(ramin, decmax, geom.degrees)
rect = [radec_NE, radec_NW, radec_SW, radec_SE]

tpatches = skymap.findTractPatchList(rect)

# region ----------------------------------------------------------------------------------------------
# iterating over tract and patches to build diaForced catalogues
if '/full_table_forced' not in diaSrc_store.keys():
    if os.path.isfile( f'{SCRATCH}/results/forced_metadata.csv' ):
        metadata = pd.read_csv(  f'{SCRATCH}/results/forced_metadata.csv')
    else:
        metas = []
        for tract, patches in tpatches:
            tract_info = skymap[tract.getId()]
            for patch in patches:
                # identify the t+p
                patch_i, patch_j = patch.getIndex()
                patch_str = '{},{}'.format(patch_i, patch_j)
                tpId = {'tract': tract.getId(), 'patch': patch_str}
                print('starting with ', tpId)
                metadata = diabutler.queryMetadata('deepDiff_forced_diaSrc', metacols,dataId=tpId)
                metadata = pd.DataFrame(metadata, columns=metacols)
                metadata = metadata[metadata['filter']!='u']
                metadata = metadata[metadata['filter']!='y']
                #metadata['tract'] = tpId['tract']
                #metadata['patch'] = tpId['patch']
                #metadata = optimize_df(metadata)
                metas.append(metadata)
        metadata = pd.concat(metas).drop_duplicates()
        metadata.to_csv( f'{SCRATCH}/results/forced_metadata.csv', index=False)

    cats = []
    path = os.path.join(forcerepo, 'deepDiff')
    diffpath = 'v{}-f{}/{}/diaForced_{}-{}-{}-{}-det{}.fits'
    for idx, idr, vn, fna, raf, detN, det in metadata.itertuples():
        pp = diffpath.format(str(vn).zfill(8), fna, raf, str(vn).zfill(8), 
                            fna, raf, detN, str(det).zfill(3))
        dpath = os.path.join(path, pp)
        if os.path.exists(dpath):
            cat_raw = diabutler.get('deepDiff_forced_diaSrc', visit=vn, 
                                    detector=det)#.asAstropy()

            if len(cat_raw) != 0:
                info = diabutler.get('calexp', visit=vn,
                                    detector=det).getInfo().getVisitInfo()
                pcal = diabutler.get('deepDiff_differenceExp_photoCalib', 
                                    visit=vn, detector=det)
                cat_cal = pcal.calibrateCatalog(cat_raw)
                catalog = cat_cal.asAstropy()
                catalog['visit_n']  = vn
                catalog['filter']   = fna
                catalog['raft']     = raf
                catalog['sensor']   = detN
                catalog['detector'] = det
                catalog['mjd']      = info.getDate().get()
                catalog['airmass']  = info.getBoresightAirmass()

                cats.append(catalog)

    mastercat = vstack(cats)
    mastercat.write(f'{SCRATCH}/results/full_table_forced_diasrc.fits')
    mastercat = mastercat.to_pandas()
    #mastercat = optimize_df(mastercat)

    fullschema = list(cat_cal.schema)
    with open('schema_forced.txt', 'w') as schemafile:
        schemafile.write('Name'.ljust(54)+' Doc\n')
        for asch in fullschema:
            field = asch.getField()
            schemafile.write(field.getName().ljust(54))
            schemafile.write(' ')
            schemafile.write(field.getDoc())
            schemafile.write('\n')

    diaSrc_store['full_table_forced'] = mastercat
    diaSrc_store.flush()
#endregion -------------------------------------------------------------------------------------------

    del(mastercat)
    del(cats)


#region  -----------------------------------------------------------------------
## iterating over every tract and patch
if 'diaObject_table' not in diaSrc_store.keys() or '/assoc_table' not in diaSrc_store.keys():
    assoc_table = []
    diaObject_table = []
    for tract, patches in tpatches:
        tract_info = skymap[tract.getId()]
        for patch in patches:
            # identify the t+p
            patch_i, patch_j = patch.getIndex()
            patch_str = '{},{}'.format(patch_i, patch_j)
            tpId = {'tract': tract.getId(), 'patch': patch_str}
            ## ask for diaObject
            try:
                dOtab = diabutler.get('deepDiff_diaObject', dataId=tpId).asAstropy()
                dOtab['tract'] = tract.getId()
                dOtab['patch'] = patch_str
                diaObject_table.append(dOtab)
                assoc_table.append(diabutler.get('deepDiff_diaObjectId', 
                                                dataId=tpId).toDataFrame())
            except:
                print(tpId, 'failed')
                continue
    assoc_table = pd.concat(assoc_table)
    diaObject_table = vstack(diaObject_table)

    diaSrc_store['diaObject_table'] = diaObject_table.to_pandas()
    diaSrc_store['assoc_table'] = assoc_table

    diaSrc_store.flush()
#endregion


## iterating over tract and patches to build diaSrc catalogues
# region ----------------------------------------------------------------------------------------------
if '/full_table' not in diaSrc_store.keys():
    if os.path.isfile( f'{SCRATCH}/results/deepdiff_metadata.csv' ):
        metadata = pd.read_csv(f'{SCRATCH}/results/deepdiff_metadata.csv')
    else:
        metas = []
        for tract, patches in tpatches:
            tract_info = skymap[tract.getId()]
            for patch in patches:
                # identify the t+p
                patch_i, patch_j = patch.getIndex()
                patch_str = '{},{}'.format(patch_i, patch_j)
                tpId = {'tract': tract.getId(), 'patch': patch_str}
                print('starting with ', tpId)
                metadata = diabutler.queryMetadata('deepDiff_diaSrc', metacols,dataId=tpId)
                metadata = pd.DataFrame(metadata, columns=metacols)
                metadata = metadata[metadata['filter']!='u']
                metadata = metadata[metadata['filter']!='y']
                #metadata['tract'] = tpId['tract']
                #metadata['patch'] = tpId['patch']      
                metas.append(metadata)
        metadata = pd.concat(metas).drop_duplicates()
        metadata.to_csv(f'{SCRATCH}/results/deepdiff_metadata.csv', index=False)

    #firstcat = None
    cats = []
    path = os.path.join(diarepo, 'deepDiff')
    diffpath = 'v{}-f{}/{}/diaSrc_{}-{}-{}-{}-det{}.fits'
    for idx, idr, vn, fna, raf, detN, det  in metadata.itertuples():
        pp = diffpath.format(str(vn).zfill(8), fna, raf, str(vn).zfill(8), 
                            fna, raf, detN, str(det).zfill(3))
        dpath = os.path.join(path, pp)
        if os.path.exists(dpath):
            print('going with V:{} Det:{}'.format(vn, det))
            cat_raw = diabutler.get('deepDiff_diaSrc', visit=vn, 
                                    detector=det)#.asAstropy()

            if len(cat_raw) != 0:
                info = diabutler.get('calexp', visit=vn, 
                                    detector=det).getInfo().getVisitInfo()
                pcal = diabutler.get('deepDiff_differenceExp_photoCalib', 
                                    visit=vn, detector=det)
                cat_cal = pcal.calibrateCatalog(cat_raw)
                catalog = cat_cal.asAstropy()

                catalog['visit_n'] = vn
                catalog['filter'] = fna
                catalog['raft'] = raf
                catalog['sensor'] = detN
                catalog['detector'] = det
                catalog['mjd'] = info.getDate().get()
                catalog['airmass']  = info.getBoresightAirmass()

                cats.append(catalog)
        last_idx = idx
        last_vn = vn

    mastercat = vstack(cats)
    try:
        mastercat.write(f'{SCRATCH}/results/full_table_diasrc.fits')
    except:
        pass
    mastercat = mastercat.to_pandas()
    #diaSrc_store[store_key] = mastercat
    #mastercat = optimize_df(mastercat)

    diaSrc_store['full_table'] = mastercat
    diaSrc_store.flush()

fullschema = list(cat_cal.schema)
with open('schema.txt', 'w') as schemafile:
    schemafile.write('Name'.ljust(54)+' Doc\n')
    for asch in fullschema:
        field = asch.getField()
        schemafile.write(field.getName().ljust(54))
        schemafile.write(' ')
        schemafile.write(field.getDoc())
        schemafile.write('\n')

print('done storing stuff in {}'.format(store))
# endregion -------------------------------------------------------------------------------------------

diaSrc_store.close()


