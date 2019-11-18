import numpy as np
import pandas as pd

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils import angularSeparation
from collections import OrderedDict as Odict

dbname = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
ObsMetaData = ObservationMetaDataGenerator(database=dbname)

def main(ramax=58, ramin=56, decmin=-32, decmax=-31, t0=59215, tm=59945):
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
    Y = X.query('d1 < 1.75 | d2 < 1.75 | d3 < 1.75 |d4 < 1.75')#.obsHistID.size
    for afilter in np.unique(Y['filter']): 
        sub = Y[Y['filter']==afilter] 
        print(afilter, sub[['airmass', 'FWHMeff', 'fiveSigmaDepth']].describe())
        print('\n')
    Y.to_csv('./catalogs+tables/visits_from_minion.csv')


if __name__=='__main__':
    import argparse
    DESC = """This script will query the minion database and search in a ra-dec
              rectangle, and time window for visits. After that it will print 
              out some statistics on their seeing size, airmass and five sigma
              depth for the group of visits that could be included in that area.
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
