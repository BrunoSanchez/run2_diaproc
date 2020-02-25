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

from glob import glob

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
assocrepo = diarepo + '/rerun/assoc_secrun'
forcerepo = assocrepo + '/rerun/forcedPhot' 
tmprepo = template_repo + '/rerun/multiband'

diabutler = Butler(forcerepo)

truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58.0_56.0_-31.0_-32.0.csv')
sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58.0_56.0_-31.0_-32.0.csv')

#region --------- we are going to clean the tables using the visits in forced repo
files = glob(diarepo+'/deepDiff/*')
visitn = []
for afile in files:
    vname = os.path.basename(afile)
    visitn.append(int(vname[1:-3]))
visit_n = np.array(visitn)

visit_filter = truth_lightc.visitn.isin(visit_n)
# tab
truth_lightc = truth_lightc[visit_filter]
sntab = sntab[sntab.N_trueobserv!=0]
#endregion ---------------------------------------------------------------------

#truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58_56_-31_-32.csv')
#sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58_56_-31_-32.csv')

#diaSrc_store = pd.HDFStore('/global/cscratch1/sd/bos0109/diaSrc_forced_fulltables_v4.h5')
diaSrc_store = pd.HDFStore('/global/homes/b/bos0109/run2_diaproc/results/diaSrc_secrun_fulltables_v4.h5')
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

# I have found out that every t+p contains the same information, a single 
# table of length 947602
#store_key = str(tract.getId())+'_'+str(patch_i)+str(patch_j)
diaSrcs_tab = diaSrc_store['full_table']
diaSrc_store['new_table'] = diaSrcs_tab

diaSrcs_tab = diaSrc_store['new_table']
diaSrcs_tab['epoch_matched'] = False
diaSrc_store.flush()
#region  -----------------------------------------------------------------------
## iterating over every tract and patch
d_tol = 2.5*u.arcsec
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
            print(tpId, 'failed \n')
            continue
assoc_table = pd.concat(assoc_table)
diaObject_table = vstack(diaObject_table) 

sn_coord = SkyCoord(ra=sntab.snra_in*u.deg, 
                    dec=sntab.sndec_in*u.deg, frame='icrs')
diaO_coord = SkyCoord(ra=diaObject_table['coord_ra'], 
                      dec=diaObject_table['coord_dec'], frame='icrs')
#endregion ---------------------------------------------------------------------
#region  -----------------------------------------------------------------------
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
sntab['matched'] = match
sntab['match_ang_dist'] = d2d.to(u.arcsec)
sntab['dia_row'] = idx
sntab['dia_id'] = diaObject_table[idx]['id']
sntab['tract'] = tpId['tract']
sntab['patch'] = tpId['patch']
sntab['n_dia_detections'] = 0

diaObject_table['match'] = matchO
diaObject_table['sn_row'] = idx_
diaObject_table['match_ang_dist'] = d2d_.to(u.arcsec)
diaObject_table['sn_id'] = sntab.iloc[idx_]['galaxy_id'].values

N_matches = np.sum(match)
#endregion ---------------------------------------------------------------------
#region  -----------------------------------------------------------------------
# making the epoch-by-epoch matching

# object from diaObject not matched would be bogus
#bogus_obj = diaObject_table[~diaObject_table['match']]
# sn not matched would be missed targets
#missed_sn = sntab[~sntab.matched]
# sn and diaObjects matched would be candidates to TP
cand_obj = diaObject_table[diaObject_table['match']]
for ic in range(len(cand_obj)):
    # both data rows
    diaC = cand_obj[ic]
    snC = sntab[sntab.dia_id==diaC['id']]
    # search for dia epochs:
    diaSrcs_ids = assoc_table[assoc_table['diaObjectId']==diaC['id']]
    # query for each epoch
    # take care of u, and y filter empty rows
    dia_lc = []
    for anid in diaSrcs_ids['diaSrcIds'].values[0]:
        srcepoch = diaSrcs_tab.query('id == {}'.format(int(anid)))
        if len(srcepoch) is not 0:
            dia_lc.append(srcepoch)
    #if len(dia_lc) is not 0:
    #    dia_lc = pd.concat(dia_lc)
    
    # search for SN epochs
    snC_lc = get_truth_LC(truth_lightc, snC.snid_in.values[0])
    sn_epoch_match = np.repeat(False, len(snC_lc))
    for a_dia_epoch in dia_lc:
    #for iepoch, a_dia_epoch in dia_lc.iterrows():
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
    sntab.loc[sntab.dia_id==diaC['id'], 'n_dia_detections'] = sn_N_detects

diaSrc_store['new_table'] = diaSrcs_tab
diaSrc_store.flush()
print(N_matches, len(sntab), N_matches/len(sntab))
#endregion ---------------------------------------------------------------------
#endregion ---------------------------------------------------------------------
diaObject_table.write('results/diaObject_table.csv', format='csv', overwrite=True)                               
sntab.to_csv('results/sntab_matched.csv')

#region --------------------------- bringing everything to same ra-dec frame ---
sq = (sntab.snra_in <= 58) & (sntab.snra_in >= 56) 
sq = sq & (sntab.sndec_in >= -32) & (sntab.sndec_in <= -31) 
sntab = sntab[sq]

diaObject_table['coord_ra_deg'] = diaObject_table['coord_ra'].to(u.deg)
diaObject_table['coord_dec_deg'] = diaObject_table['coord_dec'].to(u.deg)
sq = (diaObject_table['coord_ra_deg'] >=  56) 
sq = sq & (diaObject_table['coord_ra_deg'] <= 58)
sq = sq & (diaObject_table['coord_dec_deg'] >= -32) 
sq = sq & (diaObject_table['coord_dec_deg'] <= -31)
diaObject_table = diaObject_table[sq]

diaSrcs_tab['coord_ra_deg'] = np.rad2deg(diaSrcs_tab['coord_ra'])
diaSrcs_tab['coord_dec_deg'] = np.rad2deg(diaSrcs_tab['coord_dec'])
sq = (diaSrcs_tab['coord_ra_deg'] >=  56) 
sq = sq & (diaSrcs_tab['coord_ra_deg'] <= 58)
sq = sq & (diaSrcs_tab['coord_dec_deg'] >= -32) 
sq = sq & (diaSrcs_tab['coord_dec_deg'] <= -31)
diaSrcs_tab = diaSrcs_tab[sq]

diaSrc_store['matched_tab'] = diaSrcs_tab 
#endregion  --------------------------------------------------------------------

#region ------------------------------------ unfolding the lightcurves ---------
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
#endregion ---------------------------------------------------------------------

#region  -----------------------------------------------------------------------
## ========================================================================== ##
## ================== Plots of Matched vs NOT-matched ======================= ##
## ========================================================================== ##

#region -------------------------------------------sn tab match vs nomatch -----
ff = sntab.matched
# z
plt.figure(figsize=(12, 8))
plt.subplot(3, 3, 1)
bins = np.arange(sntab.z_in.min(), sntab.z_in.max(), 0.05)
plt.hist(sntab['z_in'], bins=bins, color='black', 
         histtype='step', label='Full')
plt.hist(sntab[ff]['z_in'], bins=bins, color='red', 
         histtype='step', label='matched')
plt.xlabel('z_in')
# mB
plt.subplot(3, 3, 2)
bins = np.arange(sntab.mB.min(), sntab.mB.max(), 0.5)
plt.hist(sntab['mB'], bins=bins, color='black', histtype='step')
plt.hist(sntab[ff]['mB'], bins=bins, color='red', histtype='step')
plt.xlabel('mB')
# c_in
plt.subplot(3, 3, 3)
bins = np.arange(sntab.c_in.min(), sntab.c_in.max(), 0.025)
plt.hist(sntab['c_in'], bins=bins, color='black', histtype='step')
plt.hist(sntab[ff]['c_in'], bins=bins, color='red', histtype='step')
plt.xlabel('c_in')
# x1_in
plt.subplot(3, 3, 4)
bins = np.arange(sntab.x1_in.min(), sntab.x1_in.max(), 0.5)
plt.hist(sntab['x1_in'], bins=bins, color='black', histtype='step')
plt.hist(sntab[ff]['x1_in'], bins=bins, color='red', histtype='step')
plt.xlabel('x1_in')
# t0_in
plt.subplot(3, 3, 5)
bins = np.arange(sntab.t0_in.min(), sntab.t0_in.max(), 100)
plt.hist(sntab['t0_in'], bins=bins, color='black', 
         histtype='step', log=False, label='Full')
plt.hist(sntab[ff]['t0_in'], bins=bins, color='red', 
         histtype='step', log=False, label='matched')
plt.xlabel('t0_in')
plt.legend(loc='best')
# coords
plt.subplot(3, 3, 6)
plt.scatter(sntab['snra_in'], sntab['sndec_in'],
            color='black', marker='x', label='Full', alpha=0.6)
plt.scatter(sntab[ff]['snra_in'], sntab[ff]['sndec_in'],
            color='red', marker='.', label='matched', alpha=0.3)
plt.xlabel('snra_in')
plt.ylabel('sndec_in')
#plt.legend(loc='best')
# x0_in
plt.subplot(3, 3, 7)
bins = np.arange(sntab.x0_in.min(), sntab.x0_in.max(), 0.5e-4)
plt.hist(sntab['x0_in'], bins=bins, color='black', 
         histtype='step', log=True, label='Full')
plt.hist(sntab[ff]['x0_in'], bins=bins, color='red', 
         histtype='step', log=True, label='matched')
plt.xlabel('x0_in')
#plt.legend(loc='best')
# mB vs z
plt.subplot(3, 3, 8)
plt.scatter(sntab['z_in'], sntab['mB'],
            color='black', marker='.', label='Full', alpha=0.6)
plt.scatter(sntab[ff]['z_in'], sntab[ff]['mB'],
            color='red', marker='x', label='matched', alpha=0.3)
plt.xlabel('z_in')
plt.ylabel('mB')
#plt.legend(loc='best')
# mB vs z
plt.subplot(3, 3, 9)
plt.scatter(sntab['mB'], sntab['t0_in'],
            color='black', marker='.', label='Full', alpha=0.6)
plt.scatter(sntab[ff]['mB'], sntab[ff]['t0_in'],
            color='red', marker='x', label='matched', alpha=0.3)
plt.xlabel('mB')
plt.ylabel('t0_in')
#plt.legend(loc='best')
plt.tight_layout()
plt.savefig('match_sntab.png', dpi=480)
plt.close()
#endregion ---------------------------------------------------------------------

#region ---------------------------------------dia Object match vs NOmatch------
ff = diaObject_table['match'].data
plt.figure(figsize=(12, 8))
plt.subplot(221)
plt.hist(diaObject_table['match'].data.astype(int), log=True)
plt.xlabel('Matched ?')
plt.subplot(222)
bins=np.logspace(0, np.log10(np.max(diaObject_table['nobs'])), num=20)
plt.hist(diaObject_table[ff]['nobs'], color='black', 
         label='matched', bins=bins, histtype='step', log=True)
plt.hist(diaObject_table[~ff]['nobs'], color='red', label='FP', 
         bins=bins, histtype='step', log=True)
plt.xlabel('N observations')
plt.xscale('log')
plt.subplot(223)
bins=np.logspace(0, np.log10(np.max(diaObject_table['match_ang_dist'])), num=20)
plt.hist(diaObject_table[ff]['match_ang_dist'], color='black', 
         label='matched', histtype='step', bins=bins, log=True)
plt.hist(diaObject_table[~ff]['match_ang_dist'], color='red', 
         label='FP', histtype='step', bins=bins, log=True)
plt.xlabel('ang dist [arcsec]')
plt.xscale('log')
plt.legend(loc='best')
plt.subplot(224)
plt.plot(np.rad2deg(diaObject_table[~ff]['coord_ra'].data), 
         np.rad2deg(diaObject_table[~ff]['coord_dec'].data), '.', color='red', 
         label='FP', alpha=0.1)
plt.plot(np.rad2deg(diaObject_table[ff]['coord_ra'].data), 
         np.rad2deg(diaObject_table[ff]['coord_dec'].data), 'x', color='black', 
         label='matched', alpha=1)
plt.vlines(x=[56., 58], ymin=-32., ymax=-31, color='black')
plt.hlines(y=[-31., -32], xmin=56., xmax=58, color='black')
plt.xlabel('ra')
plt.ylabel('dec')
plt.tight_layout()
plt.savefig('diaO_table.png')
plt.close()
#endregion ---------------------------------------------------------------------
diaSrc_store['matched_tab'] = diaSrcs_tab
diaSrc_store.flush()

diaSrc_store.close()  #############################

#region  -----------------------------get the calibration zeropoints------------
# fluxes_to_calibrate = ['base_PsfFlux_instFlux',
#                        'ip_diffim_forced_PsfFlux_instFlux',
#                        'base_CircularApertureFlux_3_0_instFlux',
#                        'base_CircularApertureFlux_4_5_instFlux',
#                        'base_CircularApertureFlux_6_0_instFlux',
#                        'base_CircularApertureFlux_9_0_instFlux',
#                        'base_CircularApertureFlux_12_0_instFlux', 
#                        'base_CircularApertureFlux_17_0_instFlux', 
#                        'base_CircularApertureFlux_25_0_instFlux']
# for aflux in fluxes_to_calibrate:
#     diaSrcs_tab[aflux+'_nJy'] = np.nan 
#     diaSrcs_tab[aflux+'_nJyErr'] = np.nan 
#     diaSrcs_tab[aflux+'_calMag'] = np.nan 
#     diaSrcs_tab[aflux+'_calMagErr'] = np.nan 
#     
# for name, srccat in diaSrcs_tab.groupby(
#     ['visit_n', 'detector']):
#     vn, det = name
# 
#     photcal = diabutler.get('deepDiff_differenceExp_photoCalib', 
#                     visit=int(vn), detector=int(det))
# 
#     for id_src, asrc in srccat.iterrows():
#         for aflux in fluxes_to_calibrate:
#             flux, err = asrc[aflux], asrc[aflux+'Err'] 
#             #print(flux, err)
#             cal = photcal.instFluxToMagnitude(flux, err)
#             diaSrcs_tab.loc[id_src, aflux+'_calMag'] = cal.value
#             diaSrcs_tab.loc[id_src, aflux+'_calMagErr'] = cal.error
#             cal = photcal.instFluxToNanojansky(flux, err)
#             diaSrcs_tab.loc[id_src, aflux+'_nJy'] = cal.value
#             diaSrcs_tab.loc[id_src, aflux+'_nJyErr'] = cal.error
# diaSrc_store['matched_tab'] = diaSrcs_tab
# diaSrc_store.flush()
#diaSrc_store.close()
#endregion ---------------------------------------------------------------------

#region  --------------------------------------- analyzing brightness of objects
# plt.suplot(2, 2, 1)
# plt.hist(diaObject_table[ff]['match_ang_dist'], color='black', 
#          label='matched', histtype='step', bins=bins, log=True)
# plt.hist(diaObject_table[~ff]['match_ang_dist'], color='red', 
#          label='FP', histtype='step', bins=bins, log=True)
# plt.xlabel('ang dist [arcsec]')
# plt.legend(loc='best')

# plt.suplot(2, 2, 2)

# plt.suplot(2, 2, 3)

# plt.suplot(2, 2, 4)

# bins=np.logspace(0, np.log10(np.max(diaObject_table['match_ang_dist'])), num=20)
# plt.hist(diaObject_table[ff]['match_ang_dist'], color='black', 
#          label='matched', histtype='step', bins=bins, log=True)
# plt.hist(diaObject_table[~ff]['match_ang_dist'], color='red', 
#          label='FP', histtype='step', bins=bins, log=True)
# plt.xlabel('ang dist [arcsec]')
# plt.legend(loc='best')



#endregion ---------------------------------------------------------------------

#region  -----------------------------------------------------------------------
# #diaSrcs_tab = diaSrc_store['new_table']
# print(len(diaSrcs_tab))
# print(np.sum(diaSrcs_tab.epoch_matched), len(diaSrcs_tab), 
#       np.sum(diaSrcs_tab.epoch_matched)/len(diaSrcs_tab))
# bogus = diaSrcs_tab[~diaSrcs_tab.epoch_matched]
# reals = diaSrcs_tab[diaSrcs_tab.epoch_matched]
#endregion  --------------------------------------------------------------------



# #endregion
# cllc = snlcs.observable & snlcs.observed 
# cllc = cllc & (snlcs.filter!='y') & (snlcs.filter!='u')
# clean_snlc = snlcs[cllc]
# print(np.sum(snlcs.epoch_DIAmatch))

