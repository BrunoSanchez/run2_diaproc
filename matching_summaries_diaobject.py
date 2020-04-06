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

"""
parsing_dia_cats

Objectives
----------

Match the variable summary catalogs to diaObject table

"""


import os
os.environ['SCRATCH']='/global/cscratch1/sd/bos0109'
SCRATCH = %env SCRATCH

import sqlite3

from glob import glob

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

from astropy import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import vstack


# load the SNe
sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58.0_56.0_-31.0_-32.0.csv')
sntab = sntab[sntab.N_trueobserv!=0]
# load the Star summary
star_summ = pd.read_csv(f'{SCRATCH}/results/star_summary_stats_rect.csv')
# load the AGN summary
agn_summ  = pd.read_csv(f'{SCRATCH}/results/agn_summ_rect.csv')

# load diasrcs and objects
store = f'{SCRATCH}/results/diaSrc_thirrun_fulltables_v1.h5'
diaSrc_store = pd.HDFStore(store, mode='r+')
diaSrc_store.open()

diaObject_table = diaSrc_store['diaObject_table']
#assoc_table = diaSrc_store['assoc_table']

# validate rectangle
sq = (sntab.snra_in <= 58) & (sntab.snra_in >= 56) 
sq = sq & (sntab.sndec_in >= -32) & (sntab.sndec_in <= -31) 
sntab = sntab[sq]

sq = (star_summ.ra <= 58) & (star_summ.ra >= 56) 
sq = sq & (star_summ.dec >= -32) & (star_summ.dec <= -31) 
star_summ = star_summ[sq]

sq = (agn_summ.ra <= 58) & (agn_summ.ra >= 56) 
sq = sq & (agn_summ.dec >= -32) & (agn_summ.dec <= -31) 
agn_summ = agn_summ[sq]
# -------------------

diaObject_table['coord_ra_deg'] = np.rad2deg(diaObject_table['coord_ra'])
diaObject_table['coord_dec_deg'] = np.rad2deg(diaObject_table['coord_dec'])
sq = (diaObject_table['coord_ra_deg'] >=  56) 
sq = sq & (diaObject_table['coord_ra_deg'] <= 58)
sq = sq & (diaObject_table['coord_dec_deg'] >= -32) 
sq = sq & (diaObject_table['coord_dec_deg'] <= -31)
diaObject_table = diaObject_table[sq]

sne = sntab[['snid_in', 'snra_in', 'sndec_in']].copy()
sne.rename(columns={'snid_in': 'id', 'snra_in': 'ra', 'sndec_in': 'dec'}, inplace=True)
stars = star_summ[['id', 'ra', 'dec']].copy()
agn = agn_summ[['galaxy_id', 'ra', 'dec']].copy()
agn.rename(columns={'galaxy_id': 'id'}, inplace=True)

truth_summ = pd.concat([sne, agn, stars])
truth_summ.to_csv(f'{SCRATCH}/results/truth_summ_id_ra_dec.csv', index=False)

truth_coord = SkyCoord(ra=truth_summ.ra.values*u.deg, 
                       dec=truth_summ.dec.values*u.deg, frame='icrs')

diaO_coord = SkyCoord(ra=diaObject_table['coord_ra'].values*u.rad, 
                     dec=diaObject_table['coord_dec'].values*u.rad, frame='icrs')
#endregion ---------------------------------------------------------------------

#region  -----------------------------------------------------------------------
#produce the actual matching
d_tol = 1.*u.arcsec
idx, d2d, d3d = truth_coord.match_to_catalog_sky(diaO_coord)
idx_, d2d_, d3d_ = diaO_coord.match_to_catalog_sky(truth_coord)

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
truth_summ['matched']        = match
truth_summ['match_ang_dist'] = d2d.to(u.arcsec)
truth_summ['dia_row']        = idx
truth_summ['dia_id']         = diaObject_table['id'].values[idx]
truth_summ['n_dia_detections']= 0

diaObject_table['match']          = matchO
diaObject_table['truth_row']         = idx_
diaObject_table['match_ang_dist'] = d2d_.to(u.arcsec)

diaObject_table['truth_id'] = truth_summ['id'].values[idx_]

N_matches = np.sum(match)
print(np.sum(match), np.sum(matchO))
#endregion

truth_summ.to_csv(f'{SCRATCH}/results/truth_summ_matched.csv')
diaSrc_store['diaObject_table_matched'] = diaObject_table
diaSrc_store.flush()

# store matches all along in different components

# match with the SNe
supernovae = pd.merge(left=sntab, right=truth_summ, left_on='snid_in', right_on='id', how='left')
supernovae.to_csv(f'{SCRATCH}/results/supernovae_summary_matched.csv', index=False)
# load the Star summary
star_matched = pd.merge(left=star_summ, right=truth_summ, left_on='id', right_on='id', how='left')
star_matched.to_csv(f'{SCRATCH}/results/star_summary_matched.csv', index=False)

# load the AGN summary
agn_matched = pd.merge(left=agn_summ, right=truth_summ, left_on='galaxy_id', right_on='id', how='left')
agn_matched.to_csv(f'{SCRATCH}/results/star_summary_matched.csv', index=False)

