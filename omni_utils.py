# 2018/03/13
import os
from scipy.io import readsav
from pathlib import Path

from hatch_python_utils import time_hist2
from geospacepy import omnireader
import pandas as pd
import numpy as np


def load_omniDB():
    dir = '/SPENCEdata/Research/database/OMNI/'
    file = 'culled_OMNI_magdata.dat'
    # omnifile = Path(dir+file)

    if os.path.exists(dir+file):
        print("Opening " + file + ' ...')
        omni = readsav.read(dir+file)
        return omni
    else:
        # doesn't exist
        print("Couldn't get OMNI IDL save file!!! Returning ...")
        return

    # try:
    #     omni_path = omnifile.resolve():
    # except FileNotFoundError:
    #     # doesn't exist
    #     print("Couldn't get OMNI IDL save file!!! Returning ...")
    #     return
    # else:
    #     # exists
    #     print("Opening " + file + ' ...')
    #     omni = idlsave.read(dir+file)
    #     return omni


def nearest(items, pivot):
    # return min(items, key=lambda x: abs(x - pivot))
    return np.argmin(abs(items - pivot))


def omni_screen(dfOMNI, requiredOMNIdict):

    import operator

    def get_truth(inp, relate, cut):
        ops = {'>': operator.gt,
               '<': operator.lt,
               '>=': operator.ge,
               '<=': operator.le,
               '=': operator.eq}
        return ops[relate](inp, cut)

    nTot = dfOMNI.shape[0]
    omni_keeps = np.ones(nTot, dtype=np.bool)

    for key, val in requiredOMNIdict.items():
        tmps = get_truth(dfOMNI[key].values, val[0], val[1])
        print("{:11s} {:s}  {:.2f} : {:d}/{:d}  ({:.2f}%)".format(key,
                                                                  val[0],
                                                                  val[1],
                                                                  np.sum(tmps),
                                                                  nTot,
                                                                  np.sum(tmps)/nTot*100.))
        omni_keeps = omni_keeps & tmps

    return omni_keeps


def omni_getter(*args,
                use_1minres=False,
                merge_highres_and_hourly_dfs=True,
                get_time_history=False,
                interp_to_input_times=False,
                interpolateArgs={'method': 'time'},
                include_AL_AU=False,
                verbose=False):
    """
    omni_getter(t_start OR time_array[,t_end])
    """

    input_is_arr = False
    if len(args) == 0:
        print("Must supply either time array as first arg, or t_start and t_end as first and second args!")
    elif len(args) == 2:
        t_start = args[0]
        t_end = args[1]
    elif len(args) == 1:
        assert (isinstance(args[0], tuple) or
                isinstance(args[0], list) or
                isinstance(args[0], np.ndarray) or
                isinstance(args[0], pd.core.series.Series) or
                isinstance(args[0], pd.core.indexes.datetimes.DatetimeIndex)), "Need array-like input!"

        input_is_arr = True

        if isinstance(args[0], pd.core.series.Series):
            if verbose:
                print("Got pd.Series")
            times = pd.DatetimeIndex(args[0])
        elif (isinstance(args[0], np.ndarray) or
              isinstance(args[0], list)):
            if verbose:
                print("Got list/array of times")
            times = pd.DatetimeIndex(args[0])
        elif isinstance(args[0], pd.core.indexes.datetimes.DatetimeIndex):
            times = args[0]

        t_start = times[0].round('1s')
        t_end = times[-1].round('1s')

    timeFmt = "%Y-%m-%d %H:%M:%S"
    print("tStart, tEnd: {:s}, {:s}".format(t_start.strftime(timeFmt),
                                            t_end.strftime(timeFmt)))

    # if interp_to_input_times:
    # totdt = t_start-t_end
    # extraTid = pd.Timedelta(int(res[0])*10,unit='min')
    nHourShift = 24
    if verbose:
        print(
            "Shift start/stop time back/fwd by {:d} hours ...".format(nHourShift))

    extraTid = pd.Timedelta(nHourShift, unit='hour')
    t_start = t_start - extraTid
    t_end = t_end + extraTid

    # t_startOR, t_end,

    # t_start = data.index.min() - datetime.timedelta(1)
    # t_end = data.index.max() + datetime.timedelta(1)

    cdf_or_txt = 'txt'
    # cdf_or_txt = 'cdf'

    if use_1minres:
        res = '1min'
    else:
        res = '5min'

    # List of poss variables in omnireader.py
    omniInt = omnireader.omni_interval(t_start,
                                       t_end,
                                       res,
                                       cdf_or_txt=cdf_or_txt)
    omniInt_1hr = omnireader.omni_interval(t_start,
                                           t_end,
                                           'hourly',
                                           cdf_or_txt=cdf_or_txt)

    epochs = omniInt['Epoch']  # time array for omni highres data
    Bx, By, Bz = omniInt['BX_GSE'], omniInt['BY_GSM'], omniInt['BZ_GSM']

    AE, SymH = omniInt['AE_INDEX'], omniInt['SYM_H']

    if include_AL_AU:
        AL, AU = omniInt['AL_INDEX'], omniInt['AU_INDEX']

    vsw, psw = omniInt['flow_speed'], omniInt['Pressure']
    borovsky_reader = omnireader.borovsky(omniInt)
    borovsky = borovsky_reader()

    newell = NewellCF_calc(vsw, Bz, By)

    epochs_1hr = omniInt_1hr['Epoch']  # datetime timestamps
    F107, Kp, DST = omniInt_1hr['F10_INDEX'], omniInt_1hr['KP'], omniInt_1hr['DST']

    SW_df = pd.DataFrame(data=np.column_stack((Bz, By, Bx, AE, SymH, vsw, psw, borovsky, newell)),
                         index=epochs,
                         columns=['Bz', 'By', 'Bx', 'AE', 'SymH', 'vsw', 'psw', 'borovsky', 'newell'])

    if include_AL_AU:
        print("Adding AL, AU indices")
        SW_df['AL'] = AL
        SW_df['AU'] = AU

    SW_df_1hr = pd.DataFrame(data=np.column_stack((F107, Kp, DST)),
                             index=epochs_1hr,
                             columns=['F107', 'Kp', 'Dst'])

    if merge_highres_and_hourly_dfs:
        print("Merging high-res and hourly OMNI data ...")

        combo = SW_df.merge(SW_df_1hr, how='outer',
                            left_index=True, right_index=True)
        combo.loc[:, 'F107'] = combo.F107.interpolate(**{'method': 'time'})
        combo.loc[:, 'Kp'] = combo.Kp.interpolate(**{'method': 'time'})
        combo.loc[:, 'Dst'] = combo.Dst.interpolate(**{'method': 'time'})

        if get_time_history:
            print("Getting time history ...")
            combo = time_hist2.time_history(combo)

        if interp_to_input_times and input_is_arr:
            print("Interping to input times ...")

            combo = combo.reindex(combo.index.union(times)).interpolate(
                **interpolateArgs).reindex(times)

        return combo
    else:
        return SW_df, SW_df_1hr


def NewellCF_calc(v, bz, by):
    # v expected in km/s
    # b's expected in nT

    NCF = np.zeros_like(v)
    NCF.fill(np.nan)
    bt = np.sqrt(by**2 + bz**2)
    bztemp = bz
    bztemp[bz == 0] = .00001
    # Calculate clock angle (theta_c = t_c)
    tc = np.arctan2(by, bztemp)

    # Ryan original--why do this ?
    # neg_tc = bt*np.cos(tc)*bz < 0
    # tc[neg_tc] = tc[neg_tc] + np.pi

    # My junkness
    # nans = np.where(np.isnan(v) | np.isnan(bz) | np.isnan(by))[0]
    # if nans.size > 0:
    #   print("Have NaNs!")
    #   mask = np.zeros_like(v,dtype=np.bool)
    #   mask[nans] = True
    #   neg_tc = bt[~mask]*np.cos(tc[~mask])*bz[~mask] < 0
    #   tc[~mask[neg_tc]] = tc[~mask[neg_tc]] + np.pi
    # else:
    #   neg_tc = bt*np.cos(tc)*bz < 0
    #   tc[neg_tc] = tc[neg_tc] + np.pi

    sintc = np.abs(np.sin(tc/2.))
    # NCF = (v**1.33333)*(sintc**2.66667)*(bt**0.66667)
    NCF = (v**(4./3.))*(sintc**(8./3.))*(bt**(2./3.))

    return NCF


def get_good_newell_inds(df):
    """
    inds = get_good_newell_inds(df)
    Get OMNI indices for which Newell function is OK
    Requires that df have columns 'By','Bz','vsw', and 'newell'
    """

    limdict = dict(Bz=10**2.64,
                   By=10**2.32,
                   vsw=10**3.2)

    OKinds = df.index >= df.index.min()
    OKinds = OKinds & (~df['newell'].isna()) & (df['newell'] > 0)

    for key in limdict.keys():
        OKinds = OKinds & (df[key] <= limdict[key])

    return OKinds
