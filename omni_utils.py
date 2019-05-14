# 2018/03/13
import os
from scipy.io import readsav
from pathlib import Path
# import time_hist2
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


def omni_getter(t_start, t_end):

    # t_start = data.index.min() - datetime.timedelta(1)
    # t_end = data.index.max() + datetime.timedelta(1)

    # cdf_or_txt = 'txt'
    cdf_or_txt = 'cdf'

    # List of poss variables in omnireader.py
    omniInt = omnireader.omni_interval(t_start,
                                       t_end,
                                       '5min',
                                       cdf_or_txt=cdf_or_txt)
    omniInt_1hr = omnireader.omni_interval(t_start,
                                           t_end,
                                           'hourly',
                                           cdf_or_txt=cdf_or_txt)

    epochs = omniInt['Epoch']  # time array for omni 5min data
    epochs = omniInt['Epoch']  # time array for omni 5min data
    Bx, By, Bz = omniInt['BX_GSE'], omniInt['BY_GSM'], omniInt['BZ_GSM']

    AE, SymH = omniInt['AE_INDEX'], omniInt['SYM_H']

    vsw, psw = omniInt['flow_speed'], omniInt['Pressure']
    borovsky_reader = omnireader.borovsky(omniInt)
    borovsky = borovsky_reader()

    newell = NewellCF_calc(vsw, Bz, By)

    epochs_1hr = omniInt_1hr['Epoch']  # datetime timestamps
    F107, Kp, DST = omniInt_1hr['F10_INDEX'], omniInt_1hr['KP'], omniInt_1hr['DST']

    SW_df = pd.DataFrame(data=np.column_stack((Bz, By, Bx, AE, SymH, vsw, psw, borovsky, newell)),
                         index=epochs,
                         columns=['Bz', 'By', 'Bx', 'AE', 'SymH', 'vsw', 'psw', 'borovsky', 'newell'])
    SW_df_1hr = pd.DataFrame(data=np.column_stack((F107, Kp, DST)),
                             index=epochs_1hr,
                             columns=['F107', 'Kp', 'DST'])

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
    NCF = (v**1.33333)*(sintc**2.66667)*(bt**0.66667)
    return NCF
