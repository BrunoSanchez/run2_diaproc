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
# NEEDS TO RUN USING desc-stack

import os
import sqlite3

import numpy as np
from scipy import stats 
import matplotlib.pyplot as plt
import pandas as pd

from astropy import time
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.table import vstack

import lsst.afw.cameraGeom
import lsst.afw.geom as afwGeom
import lsst.geom as geom

from lsst.afw.geom import makeSkyWcs
from lsst.daf.persistence import Butler
from lsst.obs.lsst.imsim import ImsimMapper

from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import angularSeparation
from lsst.sims.utils import getRotSkyPos

from lsst.sims.catUtils import supernovae
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.dust import EBVbase
from lsst.sims.photUtils.Sed import Sed

from lsst.sims.photUtils.BandpassDict import BandpassDict
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5, calcMagError_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters

from collections import OrderedDict as Odict
 
LSST_BPass = BandpassDict.loadTotalBandpassesFromFiles()
#region creating the telescope camera objects
mapper = ImsimMapper()
camera = mapper.camera
trans = np.array([detector.getTransform(lsst.afw.cameraGeom.cameraSys.PIXELS,
         lsst.afw.cameraGeom.cameraSys.FIELD_ANGLE) for detector in camera])
boxes = np.array([detector.getBBox() for detector in camera])
names = np.array([detector.getName() for detector in camera])
#endregion

# path where instCat are located is /global/cscratch1/sd/descim/Run2.2i/y3-wfd/instCat/
# path where the file outputs of each image are located is in 
#  /global/cscratch1/sd/descim/Run2.2i/y3-wfd/run/outputs/
# path where the sne catalog is at /global/cscratch1/sd/descim/Run2.2i/y3-wfd/instCat/VISITN/sne_cat_VISITN.txt.gz 
# command to check the file in terminal is gunzip -c _filepath_ 

visitn = 641066
strvisit = str(visitn).zfill(8)

instcat = '/global/cscratch1/sd/descim/Run2.2i/y3-wfd/instCat/'
outputs = '/global/cscratch1/sd/descim/Run2.2i/y3-wfd/run/outputs'
sncats = '/global/cscratch1/sd/descim/Run2.2i/y3-wfd/instCat/{}/sne_cat_{}.txt.gz'
snedb = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sne_cosmoDC2_v1.1.4_MS.db'
centroid_tmpl = '/global/cscratch1/sd/descim/Run2.2i/y3-wfd/run/outputs/00641066/centroid_{}_{}_{}.txt.gz'
centroid_header = ['SourceID', 'Flux', 'RealizedFlux', 'xPix', 'yPix', 'flags', 'GalSimType']

conn = sqlite3.connect(snedb)
c = conn.cursor()

# region get SNe magnitudes from SNCosmo
def getSNCosmomags(mjd, filt, snpars_table, snid_in='MS_9940_3541'):
    asn = snpars_table.loc[snpars_table['snid_in']==snid_in]
    sn_mod = SNObject(ra=asn.snra_in[0], dec=asn.sndec_in[0])
    sn_mod.set(z=asn.z_in[0], t0=asn.t0_in[0], x1=asn.x1_in[0],
               c=asn.c_in[0], x0=asn.x0_in[0])

    # this is probably not the thing to do
    flux = sn_mod.catsimBandFlux(mjd, LSST_BPass[filt])
    mag = sn_mod.catsimBandMag(LSST_BPass[filt], mjd, flux)
    return(flux, mag)
#endregion

#region : this is me guessing the instcat column names
instnames = ['type', 'SN_ID', 'RA', 'Dec', 'mag', 'specfile', 
         'val1', 'val2', 'val3', 'val4', 'val3', 'val4', 'point', 
         'none', 'CCM', '0.0175260442', '3.1']
#endregion 
# instance catalog for the SN in this particular visit:
sn_instcat_table = pd.read_table(sncats.format(strvisit, str(visitn)), 
                      skiprows=1, sep=' ', names=instnames)

#region set of queries to the SNCOsmo SN database                      
query_tmpl = "SELECT * FROM sne_params WHERE snid_in = "
sntables = []
for arow in sn_instcat_table.itertuples():
    query = query_tmpl + "'" + str(arow.SN_ID) + "'"
    sntab = pd.read_sql_query(query, conn)
    sntables.append(sntab)
#endregion
# Supernovae in the instance catalog present in the CosmoDC2 database: 
sntables = pd.concat(sntables)

#region querying visits database
dbname = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
#dbname = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
ObsMetaData = ObservationMetaDataGenerator(database=dbname)

res = ObsMetaData.getObservationMetaData(obsHistID=visitn)[0]
parsed = res.summary['OpsimMetaData']
ditherRa = np.rad2deg(parsed['descDitheredRA'])
ditherDec = np.rad2deg(parsed['descDitheredDec'])
ditherRot = np.rad2deg(parsed['descDitheredRotTelPos'])
rotSkyPos = getRotSkyPos(ditherRa, ditherDec, res, ditherRot)

# boresight and WCS from the values above
bsight = geom.SpherePoint(ditherRa*geom.degrees, ditherDec*geom.degrees)
orient = (90-rotSkyPos)*geom.degrees
wcs_list = np.array([makeSkyWcs(t, orient, flipX=False, boresight=bsight,
                       projection='TAN') for t in trans])

#MDJ and filter from the ObsMetaData
mjd = parsed['expMJD']
filt = parsed['filter']
#endregion
# information of the visit, like coordinates and camera rotation angles
print('MJD', 'Filter', 'ditherRA', 'ditherDec', 'ditherRot', 'rotSkyPos')
print(mjd, filt, ditherRa, ditherDec, ditherRot, rotSkyPos, '\n')

#region load every centroid file and find the SN ids 
ctrid_tables = []
for detname in names:
    ctrid_file = centroid_tmpl.format(str(visitn), detname, filt)
    ctrid_tab = pd.read_csv(ctrid_file, sep='\s+', skiprows=1, 
                            names=centroid_header, low_memory=False)
    subtab = ctrid_tab.loc[ctrid_tab.SourceID.isin(sntables.snid_in.values)].copy()
    subtab['detname'] = detname
    ctrid_tables.append(subtab)
#endregion
# the SNe in the centroid files:
ctrid_sn_tab = pd.concat(ctrid_tables)

ctrid_sn_tab['sncosmo_mag'] = -10
ctrid_sn_tab['sncosmo_flux'] = -10
ctrid_sn_tab['centroid_mag'] = -10
ctrid_sn_tab['MJD'] = mjd
ctrid_sn_tab['Filter'] = filt

for asn in ctrid_sn_tab.itertuples():
    idx = ctrid_sn_tab['SourceID']==asn.SourceID
    flux, mag = getSNCosmomags(mjd, filt, snpars_table=sntables, snid_in=asn.SourceID)
    ctrid_sn_tab.loc[idx, 'sncosmo_flux'] = flux
    ctrid_sn_tab.loc[idx, 'sncosmo_mag'] = mag

    # we actually need to check SNCosmo DB
    snpars = sntables.loc[sntables['snid_in']==asn.SourceID]
    ctrid_sn_tab.loc[idx, 'sn_z_in'] = snpars.z_in[0]
    ctrid_sn_tab.loc[idx, 'sn_Ra_in'] = snpars.snra_in[0]
    ctrid_sn_tab.loc[idx, 'sn_Dec_in'] = snpars.sndec_in[0]
    ctrid_sn_tab.loc[idx, 'sn_t0_in'] = snpars.t0_in[0]
    ctrid_sn_tab.loc[idx, 'sn_mB_in'] = snpars.mB[0]
    ctrid_sn_tab.loc[idx, 'sn_x0_in'] = snpars.x0_in[0]
    ctrid_sn_tab.loc[idx, 'sn_x1_in'] = snpars.x1_in[0]
    ctrid_sn_tab.loc[idx, 'sn_c_in'] = snpars.c_in[0]
    sn_skyp = afwGeom.SpherePoint(snpars.snra_in[0], snpars.sndec_in[0], afwGeom.degrees)
            
    detname = asn.detname
    wcs = wcs_list[int(np.where(names==detname)[0])]

    sn_pixc = wcs.skyToPixel(sn_skyp)
    xsn, ysn = sn_pixc

ctrid_sn_tab['sncosmo_flux_nmgy'] = (ctrid_sn_tab['sncosmo_flux'].values*u.mgy).to(u.nmgy).value
ctrid_sn_tab['fxratio'] = ctrid_sn_tab['Flux']/ctrid_sn_tab['sncosmo_flux_nmgy']

#ctrid_sn_tab['sncosmo_delta_flux'] = ctrid_sn_tab['mag'] - ctrid_sn_tab['sncosmo_mag']
#ctrid_sn_tab['sncosmo_delta_mag'].describe()

# not every object in the instance catalog is in the centroid files
mean_dec = np.max(sntables['sndec_in'])
plt.figure(figsize=(6, 4))
plt.grid()
plt.scatter(sntables['snra_in'], sntables['sndec_in'], 
            s=6, label='SN Instance catalog', c='black')
plt.scatter(ctrid_sn_tab['sn_Ra_in'], ctrid_sn_tab['sn_Dec_in'], 
            s=3, label='Centroid File', color='grey')
plt.gca().set_aspect(1./np.cos(mean_dec))
plt.gca().invert_xaxis()
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
plt.legend(loc='best')
plt.savefig('scatter_sn_ctroid_vs_inscat.png', dpi=480)
plt.clf()

# to find them we need to make a simple query:
lostSN = sntables.loc[~sntables.snid_in.isin(ctrid_sn_tab['SourceID'])]
bins=np.arange(0, 1.3, 0.1)
plt.figure(figsize=(6,4))
plt.hist(sntables['z_in'], label='Full instance catalog', histtype='step', color='black', bins=bins)
plt.hist(lostSN['z_in'], label='SNe Not in centroid file', histtype='stepfilled', color='grey', bins=bins)
plt.hist(lostSN['z_in'], histtype='step', color='black', bins=bins)
plt.legend(loc='upper left')
plt.xlabel('Redshift')
plt.savefig('redshifts_lostSNe.png', dpi=360)
plt.clf()

bins=np.linspace(np.min(sntables['t0_in']), np.max(sntables['t0_in']), 10)
plt.figure(figsize=(6,4))
plt.hist(sntables['t0_in'], label='Full instance catalog', histtype='step', color='black', bins=bins)
plt.hist(lostSN['t0_in'], label='SNe Not in centroid file', histtype='stepfilled', color='grey', bins=bins)
plt.hist(lostSN['t0_in'], histtype='step', color='black', bins=bins)
plt.legend(loc='upper left')
plt.xlabel('t_0 [MJD]')
plt.savefig('t_knots_lostSNe.png', dpi=360)
plt.clf()


slp, intrcp, r_val, p_val, std_err = stats.linregress(ctrid_sn_tab['Flux'], 
                                                      ctrid_sn_tab['sncosmo_flux_nmgy'])
plt.plot(ctrid_sn_tab['Flux'], ctrid_sn_tab['sncosmo_flux_nmgy'], 'ro', label='Centroid files')
plt.plot(ctrid_sn_tab['Flux'], ctrid_sn_tab['Flux']*slp + intrcp, 'k-', label='Linear fit.')
plt.title(f'Slope={slp:.3e},  intrcp={intrcp:.2e},  slope^-1={(1./slp):.1f}')
plt.xlabel('Flux from centroid file [ADU]')
plt.ylabel('Flux from SNCosmo [nMgy]')
plt.grid()
plt.legend(loc='best')
plt.savefig('flux_relationship.png', dpi=480)
plt.clf()


mm, md, sd = sigma_clipped_stats(ctrid_sn_tab['fxratio'])
plt.hist(ctrid_sn_tab['fxratio'], log=True, histtype='step', color='black', lw=1.5)
plt.vlines(x=mm, ymin=.1, ymax=300, label='Mean')
plt.vlines(mm+sd, ymin=.1, ymax=300, label='1 std', linestyles='dashed')
plt.vlines(mm-sd, ymin=.1, ymax=300, linestyles='dashed')
plt.legend(loc='best')
plt.xlabel('(Centroid file Flux)/(SNCosmo flux) [ADU/nMgy]')
plt.savefig('flux_ratios.png', dpi=480)
plt.clf()

