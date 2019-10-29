# 2019/10/29
import pandas as pd
import numpy as np


def load_CHAMP(*args,
               cal_version='2',
               all_years=False,
               add_tParm=False,
               quiet=False):
    """
    load_CHAMP(years)

    all_years   : load all available years
    """
    assert cal_version == '2', "Not have other calibrations!"

    haveYears = ['2002', '2003', '2004',
                 '2005', '2006', '2007',
                 '2008', '2009']

    if len(args) == 0:
        if all_years:
            years = haveYears
        else:
            print("Have to pick years!")
            return None
    else:
        years = args[0]

    if cal_version == '2':
        basename = 'CH-ME-2-PLPT+'

    baselocalDir = '/media/spencerh/data/CHAMP/'

    # Determine if years is list or single item
    try:
        years.sort()
    except:
        years = [years]

    # Is number?
    try:
        _ = np.exp(years[0])
        yearStr = [str(year) for year in years]
    except:
        yearStr = years

    assert all([year in haveYears for year in yearStr]
               ), "Select from available years: {:s}".format(",".join(haveYears))

    maxie = []
    for year in yearStr:
        tmpdir = baselocalDir+year+'/'
        infile = basename+year+'.pd'

        if not quiet:
            print("Loading {:s} ...".format(infile))
        maxie.append(pd.read_pickle(tmpdir+infile))

    df = pd.concat(maxie)
    if add_tParm:
        if not quiet:
            print("Adding df['tParm'] (sine of DOY) ...")
        df['tParm'] = np.sin(2*np.pi*(df['doy']-79.)/365.25)

    return df
