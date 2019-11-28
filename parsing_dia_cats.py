#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parsing_dia_cats.py
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

truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58.0_56.0_-31.0_-32.0.csv')
sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58.0_56.0_-31.0_-32.0.csv')
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
diaSrcCats = {}
for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}

        metadata = diabutler.queryMetadata('deepDiff_diaSrc', metacols, dataId=tpId)
        metadata = pd.DataFrame(metadata, columns=metacols)
        #metadata = metadata[metadata['filter']!='u']
        metadata = metadata[metadata['filter']!='y']
        cats = []
        for idx, idr, vn, fna, raf, detN, det in metadata.itertuples():
            #if fna=='y' or fna=='u': continue
            try:
                catalog = diabutler.get('deepDiff_diaSrc', visit=vn, detector=det).asAstropy()
                if len(catalog) is not 0:
                    catalog['visit_n'] = vn
                    catalog['filter'] = fna
                    catalog['raft'] = raf
                    catalog['sensor'] = detN
                    catalog['detector'] = det
                    cats.append(catalog)
            except: 
                print(tpId, vn, fna, raf, det,'failed \n')
                
        mastercat = vstack(cats)
    diaSrcCats[str(tract)+'_'+str(patch_i)+str(patch_j)] = mastercat

## iterating over every tract and patch
N_matches = 0
sn_matched_tab = []
diaO_matched_tab = []
assoc_table = []
d_tol = 2.5*u.arcsec
for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}

        ## get the t+p coordinate box corners
        tp_box = tract_info.getPatchInfo(patch.getIndex()).getOuterBBox()
        tp_pos_list = tp_box.getCorners()
        # Cast to Point2D, because pixelToSky below 
        # will refuse to work with a Point2I object.
        tp_pos_list = [afwGeom.Point2D(tp) for tp in tp_pos_list]

        wcs = tract_info.getWcs()
        corners = wcs.pixelToSky(tp_pos_list)
        corners = np.array([[c.getRa().asDegrees(), c.getDec().asDegrees()] \
                            for c in corners])
        ra, dec = corners[:, 0], corners[:, 1]
        min_ra, max_ra = np.min(ra), np.max(ra)
        min_dec, max_dec = np.min(dec), np.max(dec)

        ## filter the SN tab 
        snq = 'snra_in > {} and snra_in < {} and sndec_in > {} and sndec_in < {}' 
        SNtab = sntab.query(snq.format(min_ra, max_ra, min_dec, max_dec))
        sn_coord = SkyCoord(ra=SNtab.snra_in*u.deg, dec=SNtab.sndec_in*u.deg, 
                            frame='icrs')
        ## ask for diaObject
        try:
            diaObject_table = diabutler.get('deepDiff_diaObject', dataId=tpId).asAstropy()
            assoc_table.append(diabutler.get('deepDiff_diaObjectId', 
                                             dataId=tpId).toDataFrame())
        except:
            print(tpId, 'failed \n')
            continue
        diaO_coord = SkyCoord(ra=diaObject_table['coord_ra'], 
                              dec=diaObject_table['coord_dec'],
                              frame='icrs')

        #produce the actual matching
        idx, d2d, d3d = sn_coord.match_to_catalog_sky(diaO_coord)
        idx_, d2d_, d3d_ = diaO_coord.match_to_catalog_sky(sn_coord)

        match = np.repeat(False, len(idx))
        matchO = np.repeat(False, len(idx_))
        not_matched = []
        for i in range(len(idx)):
            if i==idx_[idx[i]] and d2d[i]<d_tol and d2d_[idx[i]]<d_tol:
                match[i] = True
                matchO[idx[i]] = True
            else:
                not_matched.append([i, idx[i], idx_[idx[i]]])
        not_matched = np.array(not_matched)
        SNtab['matched'] = match
        SNtab['match_ang_dist'] = d2d.to(u.arcsec)
        SNtab['dia_row'] = idx
        SNtab['dia_id'] = diaObject_table[idx]['id']
        SNtab['tract'] = tpId['tract']
        SNtab['patch'] = tpId['patch']
        
        diaObject_table['match'] = matchO
        diaObject_table['sn_row'] = idx_
        diaObject_table['match_ang_dist'] = d2d_.to(u.arcsec)
        diaObject_table['sn_id'] = SNtab.iloc[idx_]['galaxy_id'].values
        diaObject_table['tract'] = tpId['tract']
        diaObject_table['patch'] = tpId['patch']

        print(tract, patch, np.sum(match), np.sum(matchO))
        N_matches += np.sum(match)
        
        sn_matched_tab.append(SNtab)
        diaO_matched_tab.append(diaObject_table)

sn_matched_tab = pd.concat(sn_matched_tab)
diaO_matched_tab = vstack(diaO_matched_tab)
assoc_table = pd.concat(assoc_table)
print(N_matches, len(sntab), N_matches/len(sntab))


# plots of matched vs not matched
#region  
ff = sn_matched_tab.matched
# z
plt.figure(figsize=(12, 8))
plt.subplot(2, 3, 1)
bins = np.arange(sntab.z_in.min(), sntab.z_in.max(), 0.05)
plt.hist(sn_matched_tab[ff]['z_in'], bins=bins, color='black', 
         histtype='step', label='matched')
plt.hist(sn_matched_tab[~ff]['z_in'], bins=bins, color='red', 
         histtype='step', label='NOT matched')
plt.xlabel('z_in')
# mB
plt.subplot(2, 3, 2)
bins = np.arange(sntab.mB.min(), sntab.mB.max(), 0.5)
plt.hist(sn_matched_tab[ff]['mB'], bins=bins, color='black', histtype='step')
plt.hist(sn_matched_tab[~ff]['mB'], bins=bins, color='red', histtype='step')
plt.xlabel('mB')
# c_in
plt.subplot(2, 3, 3)
bins = np.arange(sntab.c_in.min(), sntab.c_in.max(), 0.025)
plt.hist(sn_matched_tab[ff]['c_in'], bins=bins, color='black', histtype='step')
plt.hist(sn_matched_tab[~ff]['c_in'], bins=bins, color='red', histtype='step')
plt.xlabel('c_in')
# x1_in
plt.subplot(2, 3, 4)
bins = np.arange(sntab.x1_in.min(), sntab.x1_in.max(), 0.5)
plt.hist(sn_matched_tab[ff]['x1_in'], bins=bins, color='black', histtype='step')
plt.hist(sn_matched_tab[~ff]['x1_in'], bins=bins, color='red', histtype='step')
plt.xlabel('x1_in')
# x0_in
plt.subplot(2, 3, 5)
bins = np.arange(sntab.x0_in.min(), sntab.x0_in.max(), 0.5e-4)
plt.hist(sn_matched_tab[ff]['x0_in'], bins=bins, color='black', 
         histtype='step', log=True, label='matched')
plt.hist(sn_matched_tab[~ff]['x0_in'], bins=bins, color='red', 
         histtype='step', log=True, label='NOT matched')
plt.xlabel('x0_in')
plt.legend(loc='best')
# coords
plt.subplot(2, 3, 6)
plt.scatter(sn_matched_tab[ff]['snra_in'], sn_matched_tab[ff]['sndec_in'],
            color='black', marker='x', label='matched', alpha=0.6)
plt.scatter(sn_matched_tab[~ff]['snra_in'], sn_matched_tab[~ff]['sndec_in'],
            color='red', marker='.', label='NOT matched', alpha=0.3)
plt.xlabel('snra_in')
plt.ylabel('sndec_in')
#plt.legend(loc='best')
plt.tight_layout()
plt.savefig('match_sntab.png')
plt.close()
#endregion

# object from diaObject not matched would be bogus
bogus_obj = diaO_matched_tab[diaO_matched_tab['match']]
# sn not matched would be missed targets
missed_sn = sn_matched_tab[~ff]
# sn and diaObjects matched would be candidates to TP
for ic in range(N_matches):
    # both data rows
    diaC = diaO_matched_tab[ic]
    snC = sn_matched_tab[~ff].iloc[ic]

    # search for dia epochs:
    diaSrcs_ids = assoc_table[assoc_table['diaObjectId']==diaC['id']]
    



# ---------------------------------------------------------------------------- #
# # imagine we have a tract patch
# tpId = {'tract':4431, 'patch': '1,5'}


# assoc_table = diabutler.get('deepDiff_diaObjectId', dataId=tpId).toDataFrame()  # not working yet
# #pq_path = os.path.join(assocrepo, 
# #    'deepDiff/diaObject/{tract}/{patch}/diaObjectId-{tract}-{patch}.parq')
# #assoc_table = pd.read_parquet(pq_path.format(**tpId), engine='pyarrow')
# diaObject_table = diabutler.get('deepDiff_diaObject', dataId=tpId)

# metacols = ['id', 'visit', 'filter', 'raftName', 'detectorName', 'detector']
# metadata = diabutler.queryMetadata('deepDiff_diaSrc', metacols, dataId=tpId)
# metadata = pd.DataFrame(metadata, columns=metacols)
# metadata = metadata[metadata['filter']!='u']
# metadata = metadata[metadata['filter']!='y']
# cats = []
# path = os.path.join(diarepo, 'deepDiff')
# diffpath = 'v{}-f{}/{}/diaSrc_{}-{}-{}-{}-det{}.fits'
# for idx, idr, vn, fna, raf, detN, det in metadata.itertuples():
#     if fna=='y' or fna=='u': continue
    
#     # this is the path of the table in the files
#     dpath = os.path.join(path, diffpath.format(str(vn).zfill(8), fna, raf, 
#                                                str(vn).zfill(8), fna,raf, detN, 
#                                                str(det).zfill(3)))
#     if os.path.exists(dpath):
#         catalog = diabutler.get('deepDiff_diaSrc', visit=vn, detector=det).asAstropy()
#         if len(catalog) is not 0:
#             catalog['visit_n'] = vn
#             catalog['filter'] = fna
#             catalog['raft'] = raf
#             catalog['sensor'] = detN
#             catalog['detector'] = det

#             cats.append(catalog)

# mastercat = vstack(cats)
