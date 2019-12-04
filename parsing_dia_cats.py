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

diaSrc_store = pd.HDFStore('/global/cscratch1/sd/bos0109/diaSrc_fulltables.h5')
diaSrc_store.open()
metacols = ['id', 'visit', 'filter', 'raftName', 'detectorName', 'detector']

#region ------------------------------------------------------------------------
## ========================================================================== ##
## =================== Convenient functions ================================= ##
## ========================================================================== ##
def get_truth_LC(truth_tab, snid):
    sffx = ['_observable', '_observed', '_flux', '_fluxErr', '_mag', '_magErr']
    snid = str(snid)
    colset = ['mjd', 'filter', 'visitn'] + [snid+acol for acol in sffx]
    return truth_tab[colset]
#endregion  --------------------------------------------------------------------

#region  -----------------------------------------------------------------------
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

## iterating over every tract and patch
N_matches = 0
sn_matched_tab = []
diaO_matched_tab = []
assoc_table = []
d_tol = 2.5*u.arcsec
# I have found out that every t+p contains the same information, a single 
# table of length 947602
#store_key = str(tract.getId())+'_'+str(patch_i)+str(patch_j)
diaSrcs_tab = diaSrc_store['/4640_30']
diaSrc_store['new_tab'] = diaSrcs_tab
diaSrcs_tab = diaSrc_store['new_tab']
diaSrcs_tab['epoch_matched'] = False

for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}

        ## ask for diaObject
        try:
            diaObject_table = diabutler.get('deepDiff_diaObject', 
                                            dataId=tpId).asAstropy()
            assoc_table.append(diabutler.get('deepDiff_diaObjectId', 
                                             dataId=tpId).toDataFrame())
        except:
            print(tpId, 'failed \n')
            continue
        
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
        print(min_ra, max_ra, min_dec, max_dec)
        ## filter the SN tab 
        snq = 'snra_in > {} and snra_in < {} and sndec_in > {} and sndec_in < {}' 
        SNtab = sntab.query(snq.format(min_ra, max_ra, min_dec, max_dec))
        sn_coord = SkyCoord(ra=SNtab.snra_in*u.deg, dec=SNtab.sndec_in*u.deg, 
                            frame='icrs')
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
        SNtab['n_dia_detections'] = 0
        
        diaObject_table['match'] = matchO
        diaObject_table['sn_row'] = idx_
        diaObject_table['match_ang_dist'] = d2d_.to(u.arcsec)
        diaObject_table['sn_id'] = SNtab.iloc[idx_]['galaxy_id'].values
        diaObject_table['tract'] = tpId['tract']
        diaObject_table['patch'] = tpId['patch']

        print(tract, patch, np.sum(match), np.sum(matchO))
        N_matches += np.sum(match)
        
        # making the epoch-by-epoch matching

        # object from diaObject not matched would be bogus
        #bogus_obj = diaObject_table[~diaObject_table['match']]
        # sn not matched would be missed targets
        #missed_sn = SNtab[~SNtab.matched]
        # sn and diaObjects matched would be candidates to TP
        cand_obj = diaObject_table[diaObject_table['match']]
        #cand_sn = SNtab[SNtab.matched]
        current_assoc = assoc_table[-1]
        for ic in range(len(cand_obj)):
            # both data rows
            diaC = cand_obj[ic]
            snC = SNtab[SNtab.dia_id==diaC['id']]
            # search for dia epochs:
            diaSrcs_ids = current_assoc[current_assoc['diaObjectId']==diaC['id']]
            # query for each epoch
            dia_lc = []
            for anid in diaSrcs_ids['diaSrcIds'].values[0]:
                dia_lc.append(diaSrcs_tab.query('id == {}'.format(anid)))
            # search for SN epochs
            snC_lc = get_truth_LC(truth_lightc, snC.snid_in.values[0])
            sn_epoch_match = np.repeat(False, len(snC_lc))
            for a_dia_epoch in dia_lc:
                v_match = snC_lc.visitn == int(a_dia_epoch['visit_n'])
                sn_epoch_match = sn_epoch_match | v_match
                if np.sum(v_match) == 1: # a true positive!
                    diaSrcs_tab.loc[a_dia_epoch.index[0], 'epoch_matched'] = True
                elif np.sum(v_match) == 0: pass  # a false positive
                elif np.sum(v_match) > 1: pass  # multipl matches? repeated obj?

            vep_col = snC.snid_in.values[0]+'_epoch_DIAmatch'
            if vep_col in truth_lightc.columns:
                truth_lightc[vep_col] = truth_lightc[vep_col] | sn_epoch_match
            else:
                truth_lightc[vep_col] = sn_epoch_match
            sn_N_detects = np.sum(truth_lightc[vep_col])
            SNtab.loc[SNtab.dia_id==diaC['id'], 'n_dia_detections'] = sn_N_detects

        sn_matched_tab.append(SNtab)
        diaO_matched_tab.append(diaObject_table)
        #diaSrc_store[store_key] = diaSrcs_tab
        diaSrc_store['new_tab'] = diaSrcs_tab
        diaSrc_store.flush()

sn_matched_tab = pd.concat(sn_matched_tab)
diaO_matched_tab = vstack(diaO_matched_tab)
assoc_table = pd.concat(assoc_table)
print(N_matches, len(sntab), N_matches/len(sntab))
#endregion ---------------------------------------------------------------------

#region  -----------------------------------------------------------------------
## ========================================================================== ##
## ================== Plots of Matched vs NOT-matched ======================= ##
## ========================================================================== ##
 
#region -------------------------------------------sn tab match vs nomatch -----
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
#endregion ---------------------------------------------------------------------

#region ---------------------------------------dia Object match vs NOmatch------
ff = diaO_matched_tab['match'].data
plt.subplot(221)
plt.hist(diaO_matched_tab['match'].data.astype(int), log=True)
plt.xlabel('Matched ?')
plt.subplot(222)
bins=np.logspace(0, np.log10(np.max(diaO_matched_tab['nobs'])), num=20)
plt.hist(diaO_matched_tab[ff]['nobs'], color='black', 
         label='matched', bins=bins, histtype='step', log=True)
plt.hist(diaO_matched_tab[~ff]['nobs'], color='red', label='FP', 
         bins=bins, histtype='step', log=True)
plt.xlabel('N observations')
plt.subplot(223)
bins=np.logspace(0, np.log10(np.max(diaO_matched_tab['match_ang_dist'])), num=20)
plt.hist(diaO_matched_tab[ff]['match_ang_dist'], color='black', 
         label='matched', histtype='step', bins=bins, log=True)
plt.hist(diaO_matched_tab[~ff]['match_ang_dist'], color='red', 
         label='FP', histtype='step', bins=bins, log=True)
plt.xlabel('ang dist [arcsec]')
plt.legend(loc='best')
plt.subplot(224)
plt.plot(np.rad2deg(diaO_matched_tab[~ff]['coord_ra'].data), 
         np.rad2deg(diaO_matched_tab[~ff]['coord_dec'].data), '.', color='red', 
         label='FP', alpha=0.1)
plt.plot(np.rad2deg(diaO_matched_tab[ff]['coord_ra'].data), 
         np.rad2deg(diaO_matched_tab[ff]['coord_dec'].data), 'x', color='black', 
         label='matched', alpha=0.1)
plt.vlines(x=[56., 58], ymin=-32., ymax=-31, color='black')
plt.hlines(y=[-31., -32], xmin=56., xmax=58, color='black')
plt.xlabel('ra')
plt.ylabel('dec')
plt.tight_layout()
plt.savefig('diaO_table.png')
plt.close()
#endregion ---------------------------------------------------------------------

#region  -----------------------------------------------------------------------
diaSrc_tab = diaSrc_store['new_tab']
print(len(diaSrc_tab))
print(np.sum(diaSrc_tab.epoch_matched), len(diaSrc_tab), 
      np.sum(diaSrc_tab.epoch_matched)/len(diaSrc_tab))
bogus = diaSrc_tab[~diaSrc_tab.epoch_matched]
reals = diaSrc_tab[diaSrc_tab.epoch_matched]
#endregion  --------------------------------------------------------------------

# to build the missed target samples we still need to unfold the table of
# simulated lightcurves, using the column of snid_in+_epoch_DIAmatch
# to do this we would need to iter over the 
lcs = []
for asn in sntab.itertuples():
    asnid = asn.snid_in
    asn_cols = [col for col in truth_lightc.columns if asnid+'_' in col]
    translate = {}
    for acol in asn_cols:
        translate[acol] = acol[len(asnid)+1:]
    asn_cols += ['mjd', 'filter', 'visitn']
    snlightc = truth_lightc[asn_cols].copy()
    if asnid+'_epoch_DIAmatch' not in snlightc.columns:
        snlightc['epoch_DIAmatch'] = False
    else:
        translate[asnid+'_epoch_DIAmatch'] = 'epoch_DIAmatch'
    snlightc.rename(columns=translate, inplace=True)
    snlightc['SN_id'] = asnid
    lcs.append(snlightc)
snlcs = pd.concat(lcs)
snlcs.to_csv('lightcurves/sn_matched_lcs.csv')
#endregion
cllc = snlcs.observable & snlcs.observed 
cllc = cllc & (snlcs.filter!='y') & (snlcs.filter!='u')
clean_snlc = snlcs[cllc]
print(np.sum(snlcs.epoch_DIAmatch))


assoc_table = []
diaOtables = []
for tract, patches in tpatches:
    tract_info = skymap[tract.getId()]
    for patch in patches:
        # identify the t+p
        patch_i, patch_j = patch.getIndex()
        patch_str = '{},{}'.format(patch_i, patch_j)
        tpId = {'tract': tract.getId(), 'patch': patch_str}
        try:
            diaObject_table = diabutler.get('deepDiff_diaObject', 
                                            dataId=tpId).asAstropy()
            assoc_table.append(diabutler.get('deepDiff_diaObjectId', 
                                             dataId=tpId).toDataFrame())
            diaOtables.append(diaObject_table)
        except:
            print(tpId, 'failed \n')
            continue
        #diaO_coord = SkyCoord(ra=diaObject_table['coord_ra'], 
        #                      dec=diaObject_table['coord_dec'],
        #                      frame='icrs')

