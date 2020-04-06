#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  match_variable_stars.py
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
Objective:
   take the variable objects in the DC2 simulations and locate them as 
   summary table.

   

"""

import os
os.environ['SCRATCH']='/global/cscratch1/sd/bos0109'
SCRATCH = '/global/cscratch1/sd/bos0109'

import numpy as np
import pandas as pd
import sqlite3
import gc

from astropy.table import Table, QTable
from astropy.io import fits

# Databases
#region
## jim located the databases here
jim_scratch = '/global/homes/j/jchiang8/scratch/desc/Run2.2i/stellar_variability'

# -----------
# this database is the full (DC2) stellar database, indexed with healpy
stardb = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
star_conn = sqlite3.connect(stardb)

# summary truth catalog  (DC2)
# here we have variable and not variable
star_tsumm = '/global/cscratch1/sd/descim/star_truth/star_truth_summary.db'
summ_conn = sqlite3.connect(star_tsumm)

# ------------
# the table Jim created  (JIM)
star_stats = os.path.join(jim_scratch, 'merged_star_db/star_lc_stats.db')
stats_conn = sqlite3.connect(star_stats)

# ------------
# trimmed to the DC2 region  (JIM)
star_stats = os.path.join(jim_scratch, 'star_truth_tables/star_lc_stats_trimmed.db')
stats_conn = sqlite3.connect(star_stats)

# Truth catalogs using the trimmed region  (JIM)
var_truth = os.path.join(jim_scratch, 'star_truth_tables/star_variability_truth.db')
var_truth_conn = sqlite3.connect(var_truth)
#endregion ------------

## Let's make the queries and put together the tables
#region 
query = """SELECT id, ra, dec, flux_u, flux_g, flux_r, flux_i, flux_z, flux_y 
            FROM truth_summary 
            WHERE ra <= 58 
               AND ra >=56 
               AND dec >= -32 
               AND dec<=-31 
               AND is_variable==1
      """
star_summ = pd.read_sql_query(query, summ_conn)

query = """SELECT * FROM stellar_variability"""
star_lcstats = pd.read_sql_query(query, stats_conn)

merged = pd.merge(left=star_summ, right=star_lcstats, how='inner', on='id')
del(star_lcstats)
gc.collect()

query = """SELECT * FROM stellar_variability_truth 
            WHERE obsHistID < 800000 
            AND bandpass!='u'
            AND bandpass!='y'
            """
star_variability = pd.read_sql_query(query, var_truth_conn)
lightcs = pd.merge(left=merged, 
                   right=star_variability, 
                   how='inner', on='id')
del(star_variability)
gc.collect()

# write things to disk 
#endregion
merged.to_csv(f'{SCRATCH}/results/star_summary_stats_rect.csv')
lightcs.to_csv(f'{SCRATCH}/results/star_truth_lightcs_rect.csv')

Table.from_pandas(merged).write(f'{SCRATCH}/results/star_summary_stats_rect.fits')
Table.from_pandas(lightcs).write(f'{SCRATCH}/results/star_truth_lightcs_rect.fits')

#### Let's add here the AGN too
#region AGN
agn_cosmodc2 = '/global/cscratch1/sd/jchiang8/desc/sims_GCRCatSimInterface/work/2020-02-14/agn_cosmoDC2_v1.1.4.db'
agn_conn = sqlite3.connect(agn_cosmodc2)

query = """SELECT galaxy_id, htmid_8, magNorm, redshift, M_i, ra, dec FROM agn_params 
            WHERE ra <= 58 
               AND ra >=56 
               AND dec >= -32 
               AND dec<=-31 
            """
agn_summary = pd.read_sql_query(query, agn_conn)
#endregion
agn_summary.to_csv(f'{SCRATCH}/results/agn_summ_rect.csv')
Table.from_pandas(agn_summary).write(f'{SCRATCH}/results/agn_summ_rect.fits')


