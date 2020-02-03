# 2019/10/29
import pandas as pd
import numpy as np


def load_CHAMP(*args,
               cal_version='2',
               all_years=False,
               add_tParm=False,
               useColtrane=False,
               load_fixed_MLT=True,
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

    if useColtrane:
        baselocalDir = '/Data/ift/ift_romfys1/Q1/folk/spencer/Coltrane/home/Research/database/CHAMP/'
    else:
        useExternal = False
        # baselocalDir = '/media/spencerh/data/CHAMP/'
        baselocalDir = '/SPENCEdata/Research/database/CHAMP/'

    if not quiet:
        print("baselocalDir: {:s}".format(baselocalDir))

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

        if useColtrane:
            tmpdir = baselocalDir+'/'
        else:
            if useExternal:
                tmpdir = baselocalDir+year+'/'
            else:
                tmpdir = baselocalDir+'/'

        if load_fixed_MLT:
            infile = basename+year+'_FIXEDMLT.pd'
        else:
            infile = basename+year+'.pd'
            
        if not quiet:
            print("Loading {:s} ...".format(infile))

        tmpdf = pd.read_pickle(tmpdir+infile)
        origsize = tmpdf.shape[0]
        tmpdf = tmpdf.loc[~tmpdf.index.duplicated(keep='first')]
        newsize = tmpdf.shape[0]
        print("Chucked {:d} dupe inds ...".format(origsize-newsize))
        maxie.append(tmpdf)

    df = pd.concat(maxie)
    if add_tParm:
        if not quiet:
            print("Adding df['tParm'] (sine of DOY) ...")
        df['tParm'] = np.sin(2*np.pi*(df['doy']-79.)/365.25)

    return df
