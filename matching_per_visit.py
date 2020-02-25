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

from glob import glob

from astropy import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import vstack

from lsst.afw.geom import makeSkyWcs
from lsst.daf.persistence import Butler
from lsst.obs.lsst.imsim import ImsimMapper


template_repo = '/global/cscratch1/sd/bos0109/templates_rect'
diarepo = template_repo + '/rerun/diff_rect'
assocrepo = diarepo + '/rerun/assoc_secrun'
forcerepo = assocrepo + '/rerun/forcedPhot' 
tmprepo = template_repo + '/rerun/multiband'

diabutler = Butler(forcerepo)
skymap = diabutler.get('deepCoadd_skyMap')

truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58.0_56.0_-31.0_-32.0.csv')
sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58.0_56.0_-31.0_-32.0.csv')
#truth_lightc = pd.read_csv('./lightcurves/lightcurves_cat_rect_58_56_-31_-32.csv')
#sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58_56_-31_-32.csv')

#diaSrc_store = pd.HDFStore('/global/cscratch1/sd/bos0109/diaSrc_forced_fulltables_v4.h5')
diaSrc_store = pd.HDFStore('/global/homes/b/bos0109/run2_diaproc/results/diaSrc_secrun_fulltables_v4.h5')
diaSrc_store.open()
metacols = ['id', 'visit', 'filter', 'raftName', 'detectorName', 'detector']

#region ------------------------------------------------------------------------
## =================== Convenient functions ================================= ##
def get_truth_LC(truth_tab, snid):
    sffx = ('_observable', '_observed', '_flux', '_fluxErr', '_mag', '_magErr')
    snid = str(snid)
    colset = ['mjd', 'filter', 'visitn'] + [snid+acol for acol in sffx]
    tab = truth_tab[colset].copy()
    transl = {snid+acol: acol[1:] for acol in sffx}
    tab.rename(columns=transl, inplace=True)
    tab['snid_in'] = snid
    return tab[tab.observable]
#endregion  --------------------------------------------------------------------

## ======= unfold the lightcurves =============================================
snids = [acol.strip('_observed') for acol in truth_lightc.columns 
        if '_observed' in acol]
lcs = []
for asnid in snids:
    lcs.append(get_truth_LC(truth_lightc, asnid))
lcs = pd.concat(lcs)
#======= unfold the lightcurves =============================================


diasrc_tab = diaSrc_store['full_table']

d_tol = 2.5*u.arcsec
lc_list = []
diasrc_list = []
N_matches = 0
for avisit, atab in diasrc_tab.groupby('visit_n'):
    # light curve row:
    lc = lcs.loc[lcs['visitn']==avisit].copy()
    
    snlist = sntab[sntab.snid_in.isin(lc.snid_in)]
    
    if lc.empty or snlist.empty:
        continue
    #print(avisit)
    sncoords = SkyCoord(ra=snlist.snra_in*u.deg, dec=snlist.sndec_in*u.deg)

    #srctab = atab[~atab.base_PixelFlags_flag_saturated].copy()
    #srctab = srctab[~srctab.base_PixelFlags_flag_edge].copy()
    #srctab = srctab[~srctab.ip_diffim_DipoleFit_flag_classification].copy()
    srctab = atab.copy()
    diacoords = SkyCoord(ra=srctab.coord_ra*u.rad, dec=srctab.coord_dec*u.rad)
    
    idx, d2d, d3d = sncoords.match_to_catalog_sky(diacoords)
    idx_, d2d_, d3d_ = diacoords.match_to_catalog_sky(sncoords)

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
    lc['matched'] = match
    lc['match_ang_dist'] = d2d.to(u.arcsec)
    lc['dia_row'] = idx
    lc['dia_id'] = srctab.iloc[idx]['id'].values

    srctab['cxmatch'] = matchO
    srctab['sn_row'] = idx_
    srctab['match_ang_dist'] = d2d_.to(u.arcsec)
    srctab['sn_id'] = lc.iloc[idx_]['snid_in'].values

    N_matches += np.sum(match)

    lc_list.append(lc)
    diasrc_list.append(srctab)
print(N_matches)
matched_lcs = pd.concat(lc_list)
matched_diasrc = pd.concat(diasrc_list)


bandpasses = ('g', 'r', 'i', 'z')
bins=np.arange(17, 28, 0.5)
fig, axes = plt.subplots(figsize=(6,5), nrows=2, ncols=2)

for ax, band in zip(axes.flatten(), bandpasses):
    subtab = matched_lcs.loc[matched_lcs['filter']==band]

    ax.hist(subtab.mag, histtype='step', lw=1.5, color='k', bins=bins, 
            label=f'total {band} band')

    ax.hist(subtab[subtab.matched].mag, histtype='stepfilled', lw=1.5, 
            bins=bins, color='red', label=f'matched in {band} band')
    ax.legend(loc='upper left')
plt.tight_layout()
plt.savefig('matched_lcs.png')
plt.clf()




#sn_matched = pd.read_csv('lightcurves/sn_matched_lcs.csv')
#
#files = glob(forcerepo+'/deepDiff/*')
#visitn = []
#for afile in files:
#    vname = os.path.basename(afile)
#    visitn.append(int(vname[1:-3]))
#visit_n = np.array(visitn)
#
#visit_filter = sn_matched.visitn.isin(visit_n)
#sn_matched = sn_matched[visit_filter]
#
#not_observed = sn_matched[~sn_matched.observed]
#print('Checking the not observed having no matches. #matches=', np.sum(not_observed.epoch_DIAmatch))
#subtab = sn_matched[sn_matched.observed]
#print('Checking the observed #matches=', np.sum(subtab.epoch_DIAmatch))
#not_observable = subtab[~subtab.observable]
#print('Checking the NOT observABLE #matches=', np.sum(not_observable.epoch_DIAmatch))
#observable = subtab[subtab.observable]
#print('Checking the observable #matches=', np.sum(observable.epoch_DIAmatch))


