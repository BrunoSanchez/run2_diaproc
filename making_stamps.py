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

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import dm_utilities as dmu
diaSrc_store = pd.HDFStore('/global/cscratch1/sd/bos0109/diaSrc_fulltables_v2.h5')
diaSrc_store.open()

isrc = 10
diff_visit = diaSrcs_tab.iloc[isrc]
ra, dec = diff_visit['coord_ra_deg'], diff_visit['coord_dec_deg']
diff_id = {}
# We have to convert from int64 to int to get the formatting to work right in the Gen 2 template string.
diff_id['visit'] = int(diff_visit['visit_n'])
diff_id['filter'] = diff_visit['filter']
diff_id['detector'] = int(diff_visit['detector'])

skymap = diabutler.get("%s_skyMap" % dataset_type)
coadd_id = dmu.get_coadd_id_for_ra_dec(skymap, ra, dec)
coadd_id['filter'] = diff_id['filter']

dataset_type = 'deepCoadd'
coadd_cutout = dmu.make_cutout_image(diabutler, float(ra), float(dec), 
    datasetType=dataset_type, savefits='postage_src_{}_coadd.fits'.format(isrc),
    saveplot='postage_src_{}_coadd.png'.format(isrc))

coadd_cutout = dmu.make_cutout_image2(diabutler, coadd_id, float(ra), float(dec), 
    datasetType=dataset_type, savefits='postage_src_{}_coadd.fits'.format(isrc),
    saveplot='postage_src_{}_coadd.png'.format(isrc))
science_cutout = dmu.make_cutout_image2(diabutler, diff_id, ra, dec, datasetType='calexp',
    warp_to_exposure=coadd_cutout, title='science image', 
    savefits='postage_src_{}_science.fits'.format(isrc), 
    saveplot='postage_src_{}_scienc.png'.format(isrc))

