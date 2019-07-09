# #2019/07/09
# def pandas_utils():

import numpy as np
import pandas as pd


# def resample_and_interp_df__timedate(df, interpCols,
def resample_and_interp_df__timedate(*args,
                                     resampleString='500ms',
                                     resampleFunc='mean',
                                     interpolateArgs={
                                         'method': 'time', 'limit': 30},
                                     verbose=True):
    """
    resample_and_interp_df__timedate(df[,interpCols,nidx],**kwargs)
    interpCols   : Columns to be interpolated
    nidx         : pandas date_range object (default uses floor/ceil of min/max df's indices to 1-s res)
    resampleFunc : One of 'mean', 'median', 'ffill', 'pad', [...]
    """
    nArgs = len(args)
    if nArgs > 0:
        df = args[0]

    dfCols = list(df.columns)

    # Which to interpolate?
    if nArgs > 1:
        interpCols = args[1]
    else:
        print("Interpolating all columns...")
        interpCols = list(df.columns)

    # Which to interpolate?
    if nArgs > 2:
        nidx = args[2]
    else:
        print("Making new index...")
        nidx = pd.date_range(
            df.index.min().floor('1s'), df.index.max().ceil('1s'), freq=resampleString)

    # interpCols = ['east', 'out', 'B',
    #                   'B_GEO_r', 'B_GEO_theta', 'B_GEO_phi',
    #                   'Radius', 'Latitude', 'Longitude',
    #                   'posx', 'posy', 'posz',
    #                   'bModMag', 'scspeed', 'parvelocity', 'perpvelocity']

    havers = [col for col in interpCols if (col in dfCols)]

    if len(havers) > 0:
        tmpdf = df[havers].copy()

        tmpdfnew = tmpdf.reindex(nidx).copy()

        for haver in havers:
            if haver.lower() in ['longitude', 'lng']:
                if verbose:
                    print("Unwrapping {:s} ...".format(haver))
                tmpdf.loc[:, haver] = np.unwrap(
                    np.deg2rad(tmpdf.loc[:, haver]))

            if haver.lower == 'mlt':
                if verbose:
                    print("Unwrapping {:s} ...".format(haver))
                tmpdf.loc[:, haver] = np.unwrap(
                    np.deg2rad(tmpdf.loc[:, haver]*15.))

            tmpdfnew.loc[:, haver] = tmpdf[haver].reindex(df.index.union(nidx)).interpolate(
                method='time').reindex(nidx)

            if haver.lower() in ['longitude', 'lng']:
                if verbose:
                    print("Rewrapping {:s} ...".format(haver))
                tmpdfnew.loc[:, haver] = np.rad2deg(
                    (tmpdfnew.loc[:, haver] + np.pi) % (2 * np.pi) - np.pi)

            if haver.lower() == 'mlt':
                if verbose:
                    print("Rewrapping {:s} ...".format(haver))
                tmpdfnew.loc[:, haver] = np.rad2deg(
                    (tmpdfnew.loc[:, haver] + np.pi) % (2 * np.pi))/15.

    if resampleFunc.lower() == 'mean':
        df = df.resample(resampleString).mean()[1:]
    elif resampleFunc.lower() in ['pad', 'ffill']:
        df = df.resample(resampleString).pad()[1:]

    # STARTNYE
    nHere = df.shape[0]
    if len(havers) > 0:
        # df.loc[:, havers] = tmpdfnew
        commons = df.index.intersection(tmpdfnew.index)
        for haver in havers:
            df.loc[commons, haver] = tmpdfnew.loc[commons, haver]

    return df
