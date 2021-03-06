# #2019/07/09
# def pandas_utils():

import numpy as np
import pandas as pd
from .arrays import group_consecutives

# def resample_and_interp_df__timedate(df, interpCols,


def interp_over_nans(df, checkCols,
                     max_Nsec_twixt_nans=1,
                     max_Nsec_tot=5,
                     interpolateArgs={'method': 'time', 'limit': 21}):

    max_Nsec_dt = pd.Timedelta(str(max_Nsec_twixt_nans)+'s')
    max_Nsec_dtTot = pd.Timedelta(str(max_Nsec_tot)+'s')

    naninds = (df[checkCols[0]] != df[checkCols[0]]) & False
    for col in checkCols:
        naninds = naninds | df[col].isna()
    naninds = np.where(naninds)[0]

    if naninds.size == 0:
        print("No NaNs in this df (at least not in {:s})!".format(
            ", ".join(checkCols)))
        return

    if naninds.size > 0:
        # Check if good timeresolution before and after

        nanindgroups = group_consecutives(naninds)

        print("Interping over {:d} nanInd groups ({:d} total nanInds) ...".format(
            len(nanindgroups),
            naninds.size))

        for indgrp in nanindgroups:

            nInGrp = indgrp.size
            if indgrp[0] == 0:
                startInd = 0
                stopInd = 1
            else:
                startInd = indgrp[0] - 1

            if indgrp[-1] == df.shape[0]-1:
                stopInd = df.shape[0]-1
            else:
                stopInd = indgrp[-1] + 1

            # dtBef = df.iloc[ind].name-df.iloc[startInd].name
            # dtAft = df.iloc[stopInd].name-df.iloc[ind].name
            # if (dtBef <= pd.Timedelta('1s')) and (dtAft <= pd.Timedelta('1s')):

            dtTot = df.iloc[stopInd].name-df.iloc[startInd].name
            if (dtTot/nInGrp <= max_Nsec_dt) and (dtTot <= max_Nsec_dtTot):

                # print("Good, Rickie")

                df.loc[df.iloc[startInd:stopInd+1].index,
                       :] = df.iloc[startInd:stopInd+1].interpolate(**interpolateArgs)


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
                    print("Unwrapping {:s} ...".format(haver), end=' ')
                tmpdf.loc[:, haver] = np.unwrap(
                    np.deg2rad(tmpdf.loc[:, haver]))

            if haver.lower == 'mlt':
                if verbose:
                    print("Unwrapping {:s} ...".format(haver), end=' ')
                tmpdf.loc[:, haver] = np.unwrap(
                    np.deg2rad(tmpdf.loc[:, haver]*15.))

            tmpdfnew.loc[:, haver] = tmpdf[haver].reindex(df.index.union(nidx)).interpolate(
                method='time').reindex(nidx)

            if haver.lower() in ['longitude', 'lng']:
                if verbose:
                    print("Rewrapping {:s} ...".format(haver), end=' ')
                tmpdfnew.loc[:, haver] = np.rad2deg(
                    (tmpdfnew.loc[:, haver] + np.pi) % (2 * np.pi) - np.pi)

            if haver.lower() == 'mlt':
                if verbose:
                    print("Rewrapping {:s} ...".format(haver), end=' ')
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


def nan_context_getter(df, column,
                       nPlusminus=5,
                       print_group_summary=True,
                       print_context=True):
    """
    See what things look like around NaN groups in a DataFrame
    """

    duds = pd.Series(False, index=df.index)  # Assume no duds

    try:
        # _ = (e for e in column)
        # multiColumn = True

        for col in column:
            duds = duds | df[col].duplicated()

        duds = np.where(duds)[0]

    except TypeError:
        # print my_object, 'is not iterable'
        # multiColumn = False

        duds = np.where(df[columns].duplicated())[0]

    nHere = len(df.index)

    dudGroups = group_consecutives(duds,
                                   maxDiff=1,
                                   min_streak=None,
                                   do_absDiff=False,
                                   print_summary=print_group_summary,
                                   print__maxVal=nHere)

    upBound = len(df.index)

    dudInfo = np.zeros(len(dudGroups),
                       dtype={'names': ('i', 'N',
                                        'start', 'stop',
                                        'startFrac', 'stopFrac'),
                              'formats': ('i8', 'i8',
                                          'i8', 'i8',
                                          'f8', 'f8')})

    for dudI, dudG in enumerate(dudGroups):
        minner, maxer = np.min(dudG), np.max(dudG)
        lowInd = minner-nPlusminus
        highInd = maxer+nPlusminus

        dudInfo[dudI] = (dudI, len(dudG),
                         dudG[0], dudG[-1],
                         dudG[0]/np.float64(nHere), dudG[-1]/np.float64(nHere))

        if minner == maxer:
            lines = np.arange(lowInd if lowInd >= 0 else 0,
                              highInd if highInd <= upBound else upBound)
        elif (maxer - minner) < 1000:
            lines = np.arange(lowInd if lowInd >= 0 else 0,
                              highInd if highInd <= upBound else upBound)
        # print(lines, minner, maxer)
        if print_context:
            print(df.iloc[lines])

    return np.rec.array(dudInfo)
