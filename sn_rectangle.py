import numpy as np
import pandas as pd
import sqlite3

from lsst.sims.utils import angularSeparation
from lsst.sims.catUtils import supernovae
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.dust import EBVbase

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.BandpassDict import BandpassDict
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5, calcMagError_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters

from collections import OrderedDict as Odict
from astropy import time


LSST_BandPass = BandpassDict.loadTotalBandpassesFromFiles()

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
        parsed = [Odict(obsmd.summary['OpsimMetaData']) for obsmd in res \
            if obsmd.bandpass in ("g", "r", "i", "z", "y")]
        df = pd.DataFrame(parsed)
        X = df[['obsHistID', 'filter', 'FWHMeff', 'descDitheredRA', 
                'descDitheredDec', 'airmass', 'fiveSigmaDepth', 'expMJD']].copy()
        X.descDitheredRA = np.degrees(X.descDitheredRA)
        X.descDitheredDec = np.degrees(X.descDitheredDec)
        X['d1'] = angularSeparation(ramin, decmax, X.descDitheredRA.values, X.descDitheredDec.values)
        X['d2'] = angularSeparation(ramin, decmin, X.descDitheredRA.values, X.descDitheredDec.values)
        X['d3'] = angularSeparation(ramax, decmax, X.descDitheredRA.values, X.descDitheredDec.values)
        X['d4'] = angularSeparation(ramax, decmin, X.descDitheredRA.values, X.descDitheredDec.values)
        visitab = X.query('d1 < 1.75 | d2 < 1.75 | d3 < 1.75 |d4 < 1.75')#.obsHistID.size
        del(X)
        del(df)
    
    times = visitab['expMJD']
    bands = visitab['filter']
    depths= visitab['fiveSigmaDepth']
    #colnames = ['mjd', 'filter']
    data_cols = {'mjd': times, 'filter': bands}
    for i_sn, asn in sntab.iterrows():
        sn_mod = SNObject(ra=asn.snra_in, dec=asn.sndec_in)
        sn_mod.set(z=asn.z_in, t0=asn.t0_in, x1=asn.x1_in, c=asn.c_in, x0=asn.x0_in)
        sn_flxs = [sn_mod.catsimBandFlux(mjd, LSST_BandPass[filt]) for mjd, filt in zip(times, bands)]  # done
        sn_mags = [sn_mod.catsimBandMag(LSST_BandPass[filt], mjd, flx) for mjd, filt, flx in zip(times, bands, sn_flxs)]
        sn_flxe = [sn_mod.catsimBandFluxError(mjd, LSST_BandPass[filt], m5, flx) for mjd, filt, m5, flx in zip(times, bands, depths, sn_flxs)]
        sn_mage = [sn_mod.catsimBandMagError(mjd, LSST_BandPass[filt], m5, magnitude=mag) for mjd, filt, m5, mag in zip(times, bands, depths, sn_mags)]
        data_cols[asn.snid_in+'_flux'] = sn_flxs
        data_cols[asn.snid_in+'_fluxErr'] = sn_flxe
        data_cols[asn.snid_in+'_mag'] = sn_mags
        data_cols[asn.snid_in+'_magErr'] = sn_mage
        #colnames.append(asn.snid_in)
    # dat = {}
    # for aname, adata in zip(colnames, data_cols): 
    #     dat[aname] = adata 
    lightcurves = pd.DataFrame(data_cols)
    dest_lc = './lightcurves/lightcurves_cat_rect_{}_{}_{}_{}.csv'
    lightcurves.to_csv(dest_lc.format(ramax, ramin, decmax, decmin))
    dest_snfile = './catalogs+tables/supernovae_cat_rect_{}_{}_{}_{}.csv'
    sntab.to_csv(dest_snfile.format(ramax, ramin, decmax, decmin))

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

