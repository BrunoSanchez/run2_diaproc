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
from glob import glob 

import pandas as pd
import numpy as np
import sqlite3

import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom
import lsst.geom as geom

from lsst.afw.geom import makeSkyWcs
from lsst.obs.lsst.imsim import ImsimMapper
from lsst.daf.persistence import Butler

from collections import OrderedDict as Odict
from astropy import time
import dm_utilities as dmu


## start from images up to the DM catalogs from butler/Sn table

deep_diff_path = '/global/cscratch1/sd/bos0109/templates_rect/rerun/diff_rect/deepDiff'

calexprepo = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1' 
b = Butler(calexprepo)

visit_paths = glob(deep_diff_path+'/*')
for avisit in visit_paths:
    dirname = os.path.basename(avisit)
    n_visit, filt = dirname.strip('v').split('-')
    n_visit = int(n_visit)
    filt = filt[-1]
    
    raft_list = glob(avisit+'/R*')
    for araft in raft_list:
        raft_n = int(os.path.basename(araft)[1:])
        #print(raft_n)
        # have visit number and raft number
        detlist = glob(araft+'/diaSrc_*.fits')
        for adet in detlist:
            det_n = int(os.path.basename(adet)[-8:-5])
            print(n_visit, filt, raft_n, det_n)

            im = b.get('calexp', visit=n_visit, detector=det_n)
            wcs = im.getWcs()
            
    
    
    break

