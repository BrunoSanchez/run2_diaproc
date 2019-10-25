import os
import glob
import warnings
import sqlite3
import re


import numpy as np

from astropy.table import Table

import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisplay
import lsst.afw.cameraGeom as cameraGeom

from lsst.daf.persistence import Butler

from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()

import matplotlib.pyplot as plt
import matplotlib.patches as patches


repo = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1' 
b = Butler(repo)
skymap = b.get('deepCoadd_skyMap')

# creating a rectangle of 1 sq. degree for tract/patch search
radec_NE = afwGeom.SpherePoint(58, -31, afwGeom.degrees)
radec_SE = afwGeom.SpherePoint(58, -32, afwGeom.degrees)
radec_SW = afwGeom.SpherePoint(56, -32, afwGeom.degrees)
radec_NW = afwGeom.SpherePoint(56, -31, afwGeom.degrees)
rect = [radec_NE, radec_NW, radec_SW, radec_SE]

tpatches = skymap.findTractPatchList(rect)
# check the tracts+patch, print their names
for atract in tpatches:
    print(atract[0])  # prints the number of the tract
    for patches in atract[1:]:  # next things on list are the patches
        for apatch in patches:
            print(apatch)

