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
import matplotlib.pyplot as plt
import pandas as pd

from astropy import time
from astropy import units as u
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

## this is me guessing the instcat column names
names = ['type', 'SN_ID', 'RA', 'Dec', 'mag', 'specfile', 
         'val1', 'val2', 'val3', 'val4', 'val3', 'val4', 'point', 
         'none', 'CCM', '0.0175260442', '3.1']
# instance catalog for the SN in this particular visit:
sn_cats_table = pd.read_table(sncats.format(strvisit, str(visitn)), 
                      skiprows=1, sep=' ', names=names)


# set of queries to the SNCOsmo SN database                      
query_tmpl = "SELECT * FROM sne_params WHERE snid_in = "

sntables = []
for arow in sn_cats_table.itertuples():
    query = query_tmpl + "'" + str(arow.SN_ID) + "'"
    sntab = pd.read_sql_query(query, conn)
    sntables.append(sntab)
# now we have in the sntables the 
sntables = pd.concat(sntables)

## this Y3 visit is not in the calexp repo yet. We have to find the WCS by other means
#calexp = '/global/cscratch1/sd/desc/DC2/data/Run2.2i/rerun/run2.2_calexp-v1'
#butler = Butler(calexp)

# creating the telescope camera objects
mapper = ImsimMapper()
camera = mapper.camera
trans = np.array([detector.getTransform(lsst.afw.cameraGeom.cameraSys.PIXELS,
         lsst.afw.cameraGeom.cameraSys.FIELD_ANGLE) for detector in camera])
boxes = np.array([detector.getBBox() for detector in camera])
names = np.array([detector.getName() for detector in camera])


# visits database
dbname = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
#dbname = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
ObsMetaData = ObservationMetaDataGenerator(database=dbname)

res = ObsMetaData.getObservationMetaData(obsHistID=visitn)[0]
parsed = res.summary['OpsimMetaData']
ditherRa = np.rad2deg(parsed['descDitheredRA'])
ditherDec = np.rad2deg(parsed['descDitheredDec'])
ditherRot = np.rad2deg(parsed['descDitheredRotTelPos'])
rotSkyPos = getRotSkyPos(ditherRa, ditherDec, res, ditherRot)

print(ditherRa, ditherDec, ditherRot, rotSkyPos)

bsight = geom.SpherePoint(ditherRa*geom.degrees, ditherDec*geom.degrees)
orient = (90-rotSkyPos)*geom.degrees
wcs_list = np.array([makeSkyWcs(t, orient, flipX=True, boresight=bsight,
                       projection='TAN') for t in trans])
mjd = parsed['expMJD']
filt = parsed['filter']

# load every centroid file and find the SN
ctrid_tables = []
for detname in names:
    ctrid_file = centroid_tmpl.format(str(visitn), detname, filt)
    ctrid_tab = pd.read_csv(ctrid_file, sep='\s+', skiprows=1, 
                            names=centroid_header, low_memory=False)
    subtab = ctrid_tab.loc[ctrid_tab.SourceID.isin(sntables.snid_in.values)]
    subtab['detname'] = detname
    ctrid_tables.append(subtab)
ctrid_tab = pd.concat(ctrid_tables)


sn_cats_table['sncosmo_mag'] = -10
sn_cats_table['sncosmo_flux'] = -10
sn_cats_table['centroid_mag'] = -10
for asn in sntables.itertuples():
    sn_mod = SNObject(ra=asn.snra_in, dec=asn.sndec_in)
    sn_mod.set(z=asn.z_in, t0=asn.t0_in, x1=asn.x1_in,
               c=asn.c_in, x0=asn.x0_in)

    # this is probably not the thing to do
    flux = sn_mod.catsimBandFlux(mjd, LSST_BPass[filt])
    mag = sn_mod.catsimBandMag(LSST_BPass[filt], mjd, flux)
    sn_cats_table.loc[sn_cats_table['SN_ID']==asn.snid_in, 'sncosmo_flux'] = flux
    sn_cats_table.loc[sn_cats_table['SN_ID']==asn.snid_in, 'sncosmo_mag'] = mag


    # we actually need to check the centroid files
    sn_skyp = afwGeom.SpherePoint(asn.snra_in, asn.sndec_in, afwGeom.degrees)
    contain = np.array([box.contains(afwGeom.Point2I(wcs.skyToPixel(sn_skyp))) \
                    for box, wcs in zip(boxes, wcs_list)])

    if np.sum(contain)==0:
        continue
    else:
        print('found one')
        detname = names[contain][0]
        wcs = wcs_list[contain][0]
        sn_pixc = wcs.skyToPixel(sn_skyp)
        xsn, ysn = sn_pixc
        ctrid_file = centroid_tmpl.format(str(visitn), detname, filt)

        ctrid_tab = pd.read_csv(ctrid_file, sep='\s+', skiprows=1, names=centroid_header)

        fx = ctrid_tab['xPix'] > xsn - 15. 
        fx &= ctrid_tab['xPix'] < xsn + 15.
        fy = ctrid_tab['yPix'] > ysn - 15. 
        fy &= ctrid_tab['yPix'] < ysn + 15.
        ff = fx & fy
        

        break
    




sn_cats_table['sncosmo_delta_mag'] = sn_cats_table['mag'] - sn_cats_table['sncosmo_mag']
sn_cats_table['sncosmo_delta_mag'].describe()


