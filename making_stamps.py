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

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from lsst.daf.persistence import Butler
from importlib import reload

import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom
import lsst.geom as geom
from lsst.afw.geom import makeSkyWcs
from lsst.obs.lsst.imsim import ImsimMapper

from collections import OrderedDict as Odict

import dm_utilities as dmu

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
#sntab = pd.read_csv('./catalogs+tables/supernovae_cat_rect_58_56_-31_-32.csv')
sntab = pd.read_csv('./results/sntab_matched.csv')
snlcs = pd.read_csv('lightcurves/sn_matched_lcs.csv')
visitab = pd.read_csv('./catalogs+tables/full_t_visits_from_minion.csv')

truth_lightc = truth_lightc[truth_lightc['filter']!='u']
truth_lightc = truth_lightc[truth_lightc['filter']!='y']
visitab = visitab[visitab['filter']!='u']
visitab = visitab[visitab['filter']!='y']
snlcs = snlcs[snlcs['filter']!='u']
snlcs = snlcs[snlcs['filter']!='y']

diaSrc_store = pd.HDFStore('/global/cscratch1/sd/bos0109/diaSrc_fulltables_v3.h5')
diaSrc_store.open()
diaSrcs_tab = diaSrc_store['matched_tab']
basepaths = '/global/cscratch1/sd/bos0109/run2_stamps_v3'

mapper = ImsimMapper()
camera = mapper.camera
trans = [detector.getTransform(lsst.afw.cameraGeom.cameraSys.PIXELS,
         lsst.afw.cameraGeom.cameraSys.FIELD_ANGLE) for detector in camera]
boxes = [detector.getBBox() for detector in camera]
names = [detector.getName() for detector in camera]
det_n = [detector.getId()   for detector in camera]

#region  -----------------------------------------------------------------------
stamp_path = os.path.abspath('/global/cscratch1/sd/bos0109/run2_stamps_v4/')
skymap = diabutler.get("deepCoadd_skyMap")
visit_box = Odict()
for asn in sntab[sntab.N_trueobserv>0].itertuples():
    #import ipdb; ipdb.set_trace()
    ra, dec = asn.snra_in, asn.sndec_in
    sn_skyp = afwGeom.SpherePoint(ra, dec, afwGeom.degrees)

    lc = snlcs[snlcs.SN_id==asn.snid_in]
    lc = lc[lc.observed]
    lc = lc[lc['filter']!='y']
    lc = lc[lc['filter']!='u']
    
    sndir = os.path.join(stamp_path, f'SN_stamps/{asn.snid_in}')
    if not os.path.exists(sndir):
        os.makedirs(sndir)
    head = f'id={asn.snid_in}_z={asn.z_in}_mB={asn.mB}'
    head += f'_x0={asn.x0_in}_x1={asn.x1_in}_c={asn.c_in}.snhead'
    open(os.path.join(sndir, head), 'w')

    coadd_id = dmu.get_coadd_id_for_ra_dec(skymap, ra, dec)
    for afilter, flcurve in lc.groupby('filter'):
        fpath = os.path.join(sndir, afilter)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        coadd_id['filter'] = afilter
        stamp_title = f"Coadd SN id={asn.snid_in} filter={afilter} " 
        
        coaddstamp_p = os.path.join(fpath, f'stamp_{asn.snid_in}_coadd')
        coadd_cutout = dmu.make_display_cutout_image(diabutler, coadd_id, 
            float(ra), float(dec), dataset_type='deepCoadd', title=stamp_title,
            savefits=coaddstamp_p+'.fits', saveplot=coaddstamp_p+'.png')

        for anepoch in flcurve.itertuples():
            dataId = dict(visit=anepoch.visitn)
            datarefs = list(b.subset('calexp', dataId=dataId))
            isin_some_detector=False
            # this circles through detectors
            for i, dataref in enumerate(datarefs):
                calexp = dataref.get('calexp')
                # We're not going to do anything with it here, but we can get the PSF from the calexp
                # like this:
                # psf = calexp.getPsf()
                # and we can get the zero-point (in ADU) like this
                # zero_point = calexp.getCalib().getFluxMag0()
                ccd_box = afwGeom.Box2D(calexp.getBBox())
                wcs = calexp.getWcs()
                if ccd_box.contains(wcs.skyToPixel(sn_skyp)):
                    print('sn is in dataref: ', dataref.dataId)
                    isin_some_detector=True
                    diff_id = dataref.dataId
                    break
                #center = wcs.pixelToSky(ccd_box.getCenter()).getPosition(afwGeom.degrees)
            if not isin_some_detector:
                print('no detector overlapping sn cats!') 
                continue
            #region  ----------------------------just to find detector number---
            if anepoch.visitn not in visit_box.keys():
                visitf = visitab[visitab.obsHistID==anepoch.visitn]
                if len(visitf)==0: 
                    print('visit not in table')
                    continue
                bsight = geom.SpherePoint(visitf.descDitheredRA.values[0]*geom.degrees, 
                                          visitf.descDitheredDec.values[0]*geom.degrees)
                orient = (90-visitf.descDitheredRotTelPos.values[0])*geom.degrees
    
                wcs_list = [makeSkyWcs(t, orient, flipX=False, boresight=bsight,
                                        projection='TAN') for t in trans]
                visit_box[anepoch.visitn] = [bsight, orient, wcs_list]
            else:
                bsight, orient, wcs_list = visit_box[anepoch.visitn]
            
            det_c = [(det, detn) for det, detn, box, wcs in zip(det_n, names, boxes, wcs_list) if \
                        box.contains(afwGeom.Point2I(wcs.skyToPixel(sn_skyp)))]
            if len(det_c) > 1:
                print('more than 1 detector')
                continue
            elif len(det_c)==0:
                print('no detector overlapping sn cats!') 
                continue
            else:
                detector, detname = det_c[0]
                print('detector that contains: ', detector, detname)
            import ipdb; ipdb.set_trace()
            #endregion ---------------------------------------------------------
            
            epochdir = os.path.join(fpath, f'{anepoch.visitn}')
            if not os.path.exists(epochdir):
                os.makedirs(epochdir)
            head = f'mag={anepoch.mag}_id={asn.snid_in}_z={asn.z_in}_mB={asn.mB}.epochhead'
            
            open(os.path.join(epochdir, head), 'w')                
            #diff_id = {}
            #diff_id['filter'] = afilter
            #diff_id['visit'] = int(anepoch.visitn)
            #diff_id['detector'] = int(detector)  # int(diff_visit['detector'])
            #diff_id['tract'] = coadd_id['tract']
            #diff_id['patch'] = coadd_id['patch']
            
            stamp_title = f"SN id={asn.snid_in} visit={anepoch.visitn} "
            stamp_title +=f"MJD = {anepoch.mjd} \n"
            stamp_title +=f"filter={anepoch.filter} det={detector} " 
            stamp_title +=f"matched: {anepoch.epoch_DIAmatch} " 

            scienstamp_p = os.path.join(epochdir, f'snid_{asn.snid_in}_mjd_{anepoch.mjd}_scien')
            diffstamp_p = os.path.join(epochdir, f'snid_{asn.snid_in}_mjd_{anepoch.mjd}_diff')
            try:
                science_cutout = dmu.make_display_cutout_image(b, diff_id, 
                    float(ra), float(dec), dataset_type='calexp', warp_to_exposure=coadd_cutout, 
                    title='science '+stamp_title, savefits=scienstamp_p+'.fits', 
                    saveplot=scienstamp_p+'.png')
    
                cutout_diff = dmu.make_display_cutout_image(diabutler, diff_id, 
                    float(ra), float(dec), 
                    dataset_type='deepDiff_differenceExp', warp_to_exposure=coadd_cutout, 
                    title='Diff '+stamp_title, savefits=diffstamp_p+'.fits', 
                    saveplot=diffstamp_p+'.png')
                
            except:
                print('failed the try, dataId: ', diff_id)
                continue
            
#endregion ---------------------------------------------------------------------

#region  -----------------------------------------------------------------------
for isrc in range(len(diaSrcs_tab)):
    diff_visit = diaSrcs_tab.iloc[isrc]
    ra, dec = diff_visit['coord_ra_deg'], diff_visit['coord_dec_deg']

    skymap = diabutler.get("deepCoadd_skyMap")
    coadd_id = dmu.get_coadd_id_for_ra_dec(skymap, ra, dec)
    coadd_id['filter'] = diff_visit['filter']

    diff_id = {}
    diff_id['filter'] = diff_visit['filter']
    diff_id['visit'] = int(diff_visit['visit_n'])
    diff_id['detector'] = int(diff_visit['detector'])
    diff_id['tract'] = coadd_id['tract']
    diff_id['patch'] = coadd_id['patch']
    #MJD={diff_visit['mjd']} 
    coadd_title = f"diaSrc id={diff_visit['id']} visit={diff_visit['visit_n']} "
    coadd_title +=f"filter={diff_visit['filter']} det={diff_visit['detector']} " 
    coadd_title +=f"matched: {diff_visit['epoch_matched']} " 
                   
    coaddstamp_p = basepaths + f'/stamps/diaSrc/stamp_{str(isrc).zfill(6)}_coadd'
    coadd_cutout = dmu.make_display_cutout_image(diabutler, coadd_id, 
        float(ra), float(dec), dataset_type='deepCoadd', title='Coadd '+coadd_title,
        savefits=coaddstamp_p+'.fits', saveplot=coaddstamp_p+'.png')

    scienstamp_p = basepaths + f'/stamps/diaSrc/stamp_{str(isrc).zfill(6)}_scien'
    science_cutout = dmu.make_display_cutout_image(diabutler, diff_id, 
        float(ra), float(dec), 
        dataset_type='calexp', warp_to_exposure=coadd_cutout, 
        title='science '+coadd_title, savefits=scienstamp_p+'.fits', 
        saveplot=scienstamp_p+'.png')
    
    diffstamp_p = basepaths + f'/stamps/diaSrc/stamp_{str(isrc).zfill(6)}_diff'
    cutout_diff = dmu.make_display_cutout_image(diabutler, diff_id, 
        float(ra), float(dec), 
        dataset_type='deepDiff_differenceExp', warp_to_exposure=coadd_cutout, 
        title='Diff '+coadd_title, savefits=diffstamp_p+'.fits', 
        saveplot=diffstamp_p+'.png')
#endregion ---------------------------------------------------------------------

