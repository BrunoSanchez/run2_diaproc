#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sn_rectangle.py
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

import os
import numpy as np
import pandas as pd
import sqlite3

import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom
import lsst.geom as geom

from lsst.sims.utils import angularSeparation
from lsst.sims.catUtils import supernovae
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.dust import EBVbase
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.BandpassDict import BandpassDict
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5, calcMagError_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters
from lsst.afw.geom import makeSkyWcs
from lsst.obs.lsst.imsim import ImsimMapper

from collections import OrderedDict as Odict
from astropy import time

# creating the telescope camera objects
mapper = ImsimMapper()
camera = mapper.camera
trans = [detector.getTransform(lsst.afw.cameraGeom.cameraSys.PIXELS,
         lsst.afw.cameraGeom.cameraSys.FIELD_ANGLE) for detector in camera]
boxes = [detector.getBBox() for detector in camera]
names = [detector.getName() for detector in camera]

LSST_BPass = BandpassDict.loadTotalBandpassesFromFiles()

# visits database
dbname = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
ObsMetaData = ObservationMetaDataGenerator(database=dbname)

# SNe database
database = '/global/cscratch1/sd/bos0109/sne_params_cosmoDC2_v1.1_181121.db'

conn = sqlite3.connect(database)
c = conn.cursor()

query_tmpl = "SELECT * FROM sne_params WHERE snra_in > {} AND snra_in < {} "
query_tmpl+= "AND sndec_in > {} AND sndec_in < {}"

def main(ramax=58, ramin=56, decmin=-32, decmax=-31, t0=59215, tm=61406):
    query = query_tmpl.format(ramin, ramax, decmin, decmax)

    sntab = pd.read_sql_query(query, conn)
    #sntab.to_csv('./catalogs+tables/sn_cat_rectangle.csv')

    if os.path.isfile('./catalogs+tables/full_t_visits_from_minion.csv'):
        visitab = pd.read_csv('./catalogs+tables/full_t_visits_from_minion.csv')
    else:
        res = ObsMetaData.getObservationMetaData(boundLength=2, boundType='circle', 
                                            fieldRA=(ramin-3, ramax+3), 
                                            fieldDec=(decmin-3, decmax+3), 
                                            expMJD=(t0, tm))
        parsed = [Odict(obsmd.summary['OpsimMetaData']) for obsmd in res]
        df = pd.DataFrame(parsed)
        df = df[df['filter'].isin(('g', 'r', 'i', 'z', 'y'))]

        X = df[['obsHistID', 'filter', 'FWHMeff', 'descDitheredRA', 
                'descDitheredDec', 'descDitheredRotTelPos', 'airmass', 
                'fiveSigmaDepth', 'expMJD']].copy()
        X.descDitheredRA = np.degrees(X.descDitheredRA)
        X.descDitheredDec = np.degrees(X.descDitheredDec)

        X['d1'] = angularSeparation(ramin, decmax, 
            X.descDitheredRA.values, X.descDitheredDec.values)
        X['d2'] = angularSeparation(ramin, decmin, 
            X.descDitheredRA.values, X.descDitheredDec.values)
        X['d3'] = angularSeparation(ramax, decmax, 
            X.descDitheredRA.values, X.descDitheredDec.values)
        X['d4'] = angularSeparation(ramax, decmin, 
            X.descDitheredRA.values, X.descDitheredDec.values)
        visitab = X.query('d1 < 1.75 | d2 < 1.75 | d3 < 1.75 |d4 < 1.75')
        del(X)
        del(df)
        visitab.to_csv('./catalogs+tables/full_t_visits_from_minion.csv')
    # setting the observation telescope status
    boresight = []
    orientation = []
    wcs_list = []
    for avisit in visitab.itertuples():
        bsight = geom.SpherePoint(avisit.descDitheredRA*geom.degrees, 
                                  avisit.descDitheredDec*geom.degrees)
        orient = avisit.descDitheredRotTelPos*lsst.geom.degrees
        
        wcs_list.append([makeSkyWcs(t, orient, flipX=True, boresight=bsight,
                       projection='TAN') for t in trans])
        orientation.append(orient)
        boresight.append(bsight)
    
    times = visitab['expMJD']
    bands = visitab['filter']
    depths= visitab['fiveSigmaDepth']
    #colnames = ['mjd', 'filter']
    data_cols = {'mjd': times, 'filter': bands, 'visitn': visitab['obsHistID']}
    n_observ = []
    for asn in sntab.itertuples():
        sn_mod = SNObject(ra=asn.snra_in, dec=asn.sndec_in)
        sn_mod.set(z=asn.z_in, t0=asn.t0_in, x1=asn.x1_in, 
                   c=asn.c_in, x0=asn.x0_in)

        sn_skyp = afwGeom.SpherePoint(asn.snra_in, asn.sndec_in, afwGeom.degrees)
    
        sn_flxs = []; sn_mags = []; sn_flxe = []; sn_mage = []; sn_obsrvd = []
        for mjd, filt, wcsl, m5 in zip(times, bands, wcs_list, depths):
            flux = sn_mod.catsimBandFlux(mjd, LSST_BPass[filt])
            mag = sn_mod.catsimBandMag(LSST_BPass[filt], mjd, flux)
            flux_er = sn_mod.catsimBandFluxError(mjd, LSST_BPass[filt], m5, flux)
            mag_er = sn_mod.catsimBandMagError(mjd, LSST_BPass[filt], m5, magnitude=mag)
            
            # checking sensors containing this object
            contain = [box.contains(afwGeom.Point2I(wcs.skyToPixel(sn_skyp))) \
                           for box, wcs in zip(boxes, wcsl)]
            observed = np.sum(contain) > 0
            # if observed:
            #     print('Overlaps ccd', names[np.where(contain)[0][0]])            
            sn_obsrvd.append(observed)
            sn_flxs.append(flux)  # done
            sn_mags.append(mag)
            sn_flxe.append(flux_er)
            sn_mage.append(mag_er)

        data_cols[asn.snid_in+'_observed'] = sn_obsrvd
        data_cols[asn.snid_in+'_flux'] = sn_flxs
        data_cols[asn.snid_in+'_fluxErr'] = sn_flxe
        data_cols[asn.snid_in+'_mag'] = sn_mags
        data_cols[asn.snid_in+'_magErr'] = sn_mage
        n_observ.append(np.sum(sn_obsrvd))
    sntab['Nobserv'] = n_observ
    lightcurves = pd.DataFrame(data_cols)
    dest_lc = './lightcurves/lightcurves_cat_rect_{}_{}_{}_{}.csv'
    lightcurves.to_csv(dest_lc.format(ramax, ramin, decmax, decmin))
    dest_snfile = './catalogs+tables/supernovae_cat_rect_{}_{}_{}_{}.csv'
    sntab.to_csv(dest_snfile.format(ramax, ramin, decmax, decmin))
    print("""Stored the lightcurves in {}, 
             the SN catalog in {}""".format(dest_lc.format(ramax, ramin, 
                                                           decmax, decmin), 
                                            dest_snfile.format(ramax, ramin, 
                                                               decmax, decmin)))
    return 


if __name__=='__main__':
    import argparse
    DESC = """This script will query look for supernovae by querying the 
              parameters databaes for SNe in WFD in a ra-dec rectangle. 
              After that it will also query the minion database and search for 
              visits using this rectangle and time window.
              Once this is done, it will work to figure out the specific times
              at which each SN was observed.
              """
    EPIL = """This will print a pandas describe output with columns as airmass, 
              FWHMEff and fiveSigmaDepth for each filter"""
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)

    parser.add_argument('--ramax', type=float, help='RA max [deg]')
    parser.add_argument('--ramin', type=float, help='RA min [deg]')
    parser.add_argument('--decmax', type=float, help='Dec max [deg]')
    parser.add_argument('--decmin', type=float, help='Dec min [deg]')
    parser.add_argument('--t0', type=float, help='Minimum time [MJD]', 
                        default=59215)
    parser.add_argument('--tm', type=float, help='Maximum time [MJD]',
                        default=59945)
    
    args = parser.parse_args()

    main(ramax=args.ramax, ramin=args.ramin, 
         decmin=args.decmin, decmax=args.decmax, 
         t0=args.t0, tm=args.tm)


    # plt.scatter(sntab.snra_in, sntab.sndec_in, s=20/sntab.z_in, c=sntab.z_in, cmap='RdBu_r')
    # clb = plt.colorbar()
    # clb.set_label(r'$z$', labelpad=-30, y=-0.05, rotation=0)
    # plt.xlabel(r'$\alpha$')
    # plt.ylabel(r'$\delta$')
    # t0 = time.Time(sntab.t0_in, format='mjd', scale='utc')
    # plt.figure(figsize=(12,4))
    # plt.subplot(121)
    # plt.hist(t0.datetime64, cumulative=False)
    # plt.xlabel('Year')
    # plt.ylabel('SNe peaked')

    # plt.subplot(122)
    # plt.hist(sntab.z_in, cumulative=False)
    # plt.xlabel(r'Redshift $z$')
    # plt.show()

    # plt.figure(figsize=(12,4))
    # plt.subplot(121)
    # plt.hist(sntab.x1_in, cumulative=False)
    # plt.ylabel('N')
    # plt.xlabel('SNe stretch')

    # plt.subplot(122)
    # plt.hist(sntab.c_in, cumulative=False)
    # plt.xlabel(r'Color')
    # plt.show()
