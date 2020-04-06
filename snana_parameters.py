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
os.environ['SCRATCH']='/global/cscratch1/sd/bos0109'
SCRATCH = '/global/cscratch1/sd/bos0109'

import sqlite3

import numpy as np

#import lsst.afw.cameraGeom
#import lsst.afw.geom as afwGeom   ## deprecated
import lsst.geom as geom

import matplotlib.pyplot as plt
import pandas as pd

from astropy import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import vstack
from astropy.stats import sigma_clipped_stats

from lsst.afw.geom import makeSkyWcs
from lsst.daf.persistence import Butler
from lsst.obs.lsst.imsim import ImsimMapper

from joblib import Parallel, delayed
from time import time


template_repo = '/global/cscratch1/sd/bos0109/templates_rect'
tmprepo = template_repo + '/rerun/multiband'

diarepo = template_repo + '/rerun/diff_rect'
assocrepo = diarepo + '/rerun/assoc_thirrun'
forcerepo = assocrepo + '/rerun/forcedPhot' 

diabutler = Butler(forcerepo)


if os.path.isfile( f'{SCRATCH}/results/deepdiff_metadata.csv' ):
    metadata = pd.read_csv(f'{SCRATCH}/results/deepdiff_metadata.csv')
else:
    metas = []
    for tract, patches in tpatches:
        tract_info = skymap[tract.getId()]
        for patch in patches:
            # identify the t+p
            patch_i, patch_j = patch.getIndex()
            patch_str = '{},{}'.format(patch_i, patch_j)
            tpId = {'tract': tract.getId(), 'patch': patch_str}
            print('starting with ', tpId)
            metadata = diabutler.queryMetadata('deepDiff_diaSrc', metacols,dataId=tpId)
            metadata = pd.DataFrame(metadata, columns=metacols)
            metadata = metadata[metadata['filter']!='u']
            metadata = metadata[metadata['filter']!='y']
            #metadata['tract'] = tpId['tract']
            #metadata['patch'] = tpId['patch']      
            metas.append(metadata)
    metadata = pd.concat(metas).drop_duplicates()
    metadata.to_csv(f'{SCRATCH}/results/deepdiff_metadata.csv', index=False)

rafac = 2*np.sqrt(2*np.log(2))
cats = []
path = os.path.join(diarepo, 'deepDiff')
diffpath = 'v{}-f{}/{}/diaSrc_{}-{}-{}-{}-det{}.fits'

def get_metadata(metadf):
    metacopy = metadf.copy()
    diabutler = Butler(forcerepo)
    #firstcat = None\
    for idx, idr, vn, fna, raf, detN, det  in metadf.itertuples():
        pp = diffpath.format(str(vn).zfill(8), fna, raf, str(vn).zfill(8), 
                            fna, raf, detN, str(det).zfill(3))
        dpath = os.path.join(path, pp)
        if not os.path.exists(dpath):
            continue

        metacopy.loc[idx, 'DATA_FOUND'] = True
        print('going with V:{} Det:{}'.format(vn, det))

        cal = diabutler.get('calexp', visit=vn, detector=det)
        pcal = diabutler.get('calexp_photoCalib', visit=vn, detector=det)
        cal_wcs = cal.getWcs()
        cal_box = cal.getBBox()
        cal_info = cal.getInfo().getVisitInfo()
        cal_psf = cal.getPsf()
        cal_var = cal.getVariance()
        
        mjd = cal_info.getDate().get()
        airmass = cal_info.getBoresightAirmass()

        psf_shape = cal_psf.computeShape()
        ixx, iyy, ixy = psf_shape.getParameterVector()
        # pow(ixx*iyy-ixy^2,0.25)
        psf_detradius = psf_shape.getDeterminantRadius()
        #sqrt(0.5*(ixx+iyy))
        psf_traradius = psf_shape.getTraceRadius()
        # to FWHM by multiplying by 2*math.sqrt(2*math.log(2))
        psf_detfwhm = psf_detradius*rafac
        psf_trafwhm = psf_traradius*rafac

        c1 = cal_wcs.pixelToSky(x=cal_box.beginX, y=cal_box.beginY)
        c2 = cal_wcs.pixelToSky(x=cal_box.beginX, y=cal_box.endY)
        c3 = cal_wcs.pixelToSky(x=cal_box.endX, y=cal_box.beginY)
        c4 = cal_wcs.pixelToSky(x=cal_box.endX, y=cal_box.endY)
        ra1  = c1.getRa().asDegrees()
        dec1 = c1.getDec().asDegrees()
        ra2  = c2.getRa().asDegrees()
        dec2 = c2.getDec().asDegrees()
        ra3  = c3.getRa().asDegrees()
        dec3 = c3.getDec().asDegrees()
        ra4  = c4.getRa().asDegrees()
        dec4 = c4.getDec().asDegrees()

        varPlane = cal_var.array
        sigmaPlane = np.sqrt(varPlane)
        mean_var, median_var, std_var = sigma_clipped_stats(varPlane)
        mean_sig, median_sig, std_sig = sigma_clipped_stats(sigmaPlane)

        zeroflux = pcal.getInstFluxAtZeroMagnitude()
        trsf_zflux = 2.5*np.log10(zeroflux)
        zeroflux_njy = pcal.instFluxToNanojansky(zeroflux)
        trsf_zflux_njy = 2.5*np.log10(zeroflux_njy)
        calib_mean = pcal.getCalibrationMean()
        calib_err = pcal.getCalibrationErr()
        twenty_flux = pcal.magnitudeToInstFlux(20.)

        metacopy.loc[idx, 'MJD'] = mjd
        metacopy.loc[idx, 'Airmass'] = airmass
        metacopy.loc[idx, 'PSF_Ixx'] = ixx
        metacopy.loc[idx, 'PSF_Iyy'] = iyy
        metacopy.loc[idx, 'PSF_Ixy'] = ixy 
        metacopy.loc[idx, 'PSF_detRadius'] = psf_detradius
        metacopy.loc[idx, 'PSF_traRadius'] = psf_traradius
        metacopy.loc[idx, 'PSF_detFWHM'] = psf_detfwhm
        metacopy.loc[idx, 'PSF_traFWHM'] = psf_trafwhm
        metacopy.loc[idx, 'CCD_corner_1_ra']  =  ra1
        metacopy.loc[idx, 'CCD_corner_1_dec'] = dec1
        metacopy.loc[idx, 'CCD_corner_2_ra']  =  ra2
        metacopy.loc[idx, 'CCD_corner_2_dec'] = dec2
        metacopy.loc[idx, 'CCD_corner_3_ra']  =  ra3
        metacopy.loc[idx, 'CCD_corner_3_dec'] = dec3
        metacopy.loc[idx, 'CCD_corner_4_ra']  =  ra4
        metacopy.loc[idx, 'CCD_corner_4_dec'] = dec4
        metacopy.loc[idx, 'mean_variance'] = mean_var
        metacopy.loc[idx, 'median_variance'] = median_var
        metacopy.loc[idx, 'std_variance'] = std_var
        metacopy.loc[idx, 'mean_sig'] = mean_sig
        metacopy.loc[idx, 'median_sig'] = median_sig
        metacopy.loc[idx, 'std_sig'] = std_sig
        metacopy.loc[idx, 'zeroflux'] = zeroflux
        metacopy.loc[idx, 'trsf_zflux'] = trsf_zflux
        metacopy.loc[idx, 'zeroflux_njy'] = zeroflux_njy
        metacopy.loc[idx, 'trsf_zflux_njy'] = trsf_zflux_njy
        metacopy.loc[idx, 'calib_mean'] = calib_mean
        metacopy.loc[idx, 'calib_err'] = calib_err
        metacopy.loc[idx, 'twenty_flux'] = twenty_flux
    return(metacopy)

n_jobs = 31
size = len(metadata)//n_jobs
dfs = [metadata.iloc[size*i:size*(i+1)] for i in range(n_jobs)]

print('starting the parallelization')
t0_p = time()
with Parallel(n_jobs=n_jobs, prefer='processes') as jobs:
    print('shooting jobs')
    batch_res = jobs(delayed(get_metadata)(df) for df in dfs)
print('Done calculating, concatenating')
res_p = pd.concat(batch_res)
dt_p = time() - t0_p

print('lasted ', dt_p, ' seconds')


metacopy_valid = res_p[res_p.DATA_FOUND==True]
metacopy_valid.to_csv(f'{SCRATCH}/results/ccds_metadata.csv', index=False)
print('stored stuff')