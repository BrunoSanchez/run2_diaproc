import os
import numpy as np

from astropy.table import Table

import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisplay

from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()

import matplotlib.pyplot as plt


def get_coadd_id_for_ra_dec(skymap, ra, dec):
    """
    Return a dict suitable for use as a data ID for a DM Butler

    Parameters
    ----------
    skymap: lsst.afw.skyMap.SkyMap [optional] 
        Pass in to avoid the Butler read.  Useful if you have lots of such reads.
        The skymap is just used to get the appropriate tract, patch.
        If you want to warp to a different frame, see `wcs`.
    ra: float
        Right ascension of the center of the cutout, degrees
    dec: float
        Declination of the center of the cutout, degrees
    
    Returns
    -------
    dict  - Suitable for use as a DM Butler data ID.
    """   
    # Look up the tract, patch for the RA, Dec
    radec = afwGeom.SpherePoint(ra, dec, afwGeom.degrees)
    tract_info = skymap.findTract(radec)
    patch_info = tract_info.findPatch(radec)
    coadd_id = {'tract': tract_info.getId(), 'patch': "%d,%d" % patch_info.getIndex()}

    return coadd_id


def cutout_coadd_ra_dec(butler, ra, dec, filter='r', datasetType='deepCoadd', **kwargs):
    """
    Produce a cutout from coadd from the given butler at the given RA, Dec in decimal degrees.
    
    Notes
    -----
    Trivial wrapper around 'cutout_coadd_spherepoint'
    
    Parameters
    ----------
    butler: lsst.daf.persistence.Butler
        Servant providing access to a data repository
    ra: float
        Right ascension of the center of the cutout, degrees
    dec: float
        Declination of the center of the cutout, degrees
    filter: string
        Filter of the image to load
        
    Returns
    -------
    MaskedImage
    """
    radec = afwGeom.SpherePoint(ra, dec, afwGeom.degrees)
    return cutout_coadd_spherepoint(butler, radec, filter=filter, datasetType=datasetType)
    

def cutout_coadd_spherepoint(butler, radec, filter='r', datasetType='deepCoadd',
                                  skymap=None, cutoutSideLength=51, **kwargs):
    """
    Produce a cutout from a coadd at the given afw SpherePoint radec position.
    
    Parameters
    ----------
    butler: lsst.daf.persistence.Butler
        Servant providing access to a data repository
    radec: lsst.afw.geom.SpherePoint 
        Coordinates of the center of the cutout.
    filter: string 
        Filter of the image to load
    datasetType: string ['deepCoadd']  
        Which type of coadd to load.  Doesn't support 'calexp'
    skymap: lsst.afw.skyMap.SkyMap [optional] 
        Pass in to avoid the Butler read.  Useful if you have lots of them.
    cutoutSideLength: float [optional] 
        Side of the cutout region in pixels.
    
    Returns
    -------
    MaskedImage
    """
    cutoutSize = afwGeom.ExtentI(cutoutSideLength, cutoutSideLength)

    if skymap is None:
        skymap = butler.get("%s_skyMap" % datasetType)
    
    # Look up the tract, patch for the RA, Dec
    tractInfo = skymap.findTract(radec)
    patchInfo = tractInfo.findPatch(radec)
    xy = afwGeom.PointI(tractInfo.getWcs().skyToPixel(radec))
    bbox = afwGeom.BoxI(xy - cutoutSize//2, cutoutSize)

    coaddId = {'tract': tractInfo.getId(), 'patch': "%d,%d" % patchInfo.getIndex(), 'filter': filter}
    
    cutout_image = butler.get(datasetType+'_sub', bbox=bbox, immediate=True, dataId=coaddId)
    
    return cutout_image
    
def make_cutout_image(butler, ra, dec, filter='r', vmin=None, vmax=None, label=None,
                      show=True, saveplot=False, savefits=False,
                      datasetType='deepCoadd'):
    """
    Generate and optionally display and save a postage stamp for a given RA, Dec.
    
    Parameters
    ----------
    butler: lsst.daf.persistence.Butler
        Servant providing access to a data repository
    ra: float
        Right ascension of the center of the cutout, degrees
    dec: float
        Declination of the center of the cutout, degrees
    filter: string 
        Filter of the image to load
    Returns
    -------
    MaskedImage

    Notes
    -----
    Uses matplotlib to generate stamps.  Saves FITS file if requested.
    """

    cutout_image = cutout_coadd_ra_dec(butler, ra, dec, filter=filter, datasetType='deepCoadd')
    if savefits:
        if isinstance(savefits, str):
            filename = savefits
        else:
            filename = 'postage-stamp.fits'
        cutout_image.writeFits(filename)
    
    radec = afwGeom.SpherePoint(ra, dec, afwGeom.degrees)
    xy = cutout_image.getWcs().skyToPixel(radec)
    
    if vmin is None or vmax is None:
        vmin, vmax = zscale.get_limits(cutout_image.image.array)

    plt.imshow(cutout_image.image.array, vmin=vmin, vmax=vmax, cmap='binary_r', origin='lower')
    plt.colorbar()
    plt.scatter(xy.getX() - cutout_image.getX0(), xy.getY() - cutout_image.getY0(),
                color='none', edgecolor='red', marker='o', s=200)
    if label is not None:
        plt.title(label)
    if saveplot:
        if isinstance(saveplot, str):
            filename = saveplot
        else:
            filename = 'postage-stamp.png'
        plt.savefig(filename)
    if show:
        plt.show()

    return cutout_image
    
    
def display_cutout_image(butler, ra, dec, vmin=None, vmax=None, label=None,
                      frame=None, display=None, backend='matplotlib',
                      show=True, saveplot=False, savefits=False,
                      old_matplotlib = False, filter='r',
                      datasetType='deepCoadd'):
    """
    Display a postage stamp for a given RA, Dec using LSST lsst.afw.display.
    
    Parameters
    ----------
    backend: string
        Backend can be anything that lsst.afw.display and your configuration supports: 
        e.g. matplotlib, ds9, ginga, firefly.
    
    Returns
    -------
    MaskedImage
    
    Notes
    -----
    Parameters are the same as for make_cutout_image, except for the backend.
    You definitely have the matplotlib backend.
    ds9, ginga, and firefly can be set up but are non-trivial on the scale of a simple Notebook.
    """
    cutout_image = cutout_coadd_ra_dec(butler, ra, dec, filter=filter, datasetType='deepCoadd')
    if savefits:
        if isinstance(savefits, str):
            filename = savefits
        else:
            filename = 'postage-stamp.fits'
        cutout_image.writeFits(filename)
    
    if display is None:
        display = afwDisplay.Display(frame=frame, backend=backend)

    radec = afwGeom.SpherePoint(ra, dec, afwGeom.degrees)
    xy = cutout_image.getWcs().skyToPixel(radec)
    
    display.mtv(cutout_image)
    display.scale("asinh", "zscale")
    display.dot('o', xy.getX(), xy.getY(), ctype='red')
    display.show_colorbar()
    
    return cutout_image