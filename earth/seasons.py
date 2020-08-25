import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import os


def get_season_inds(df,seasonStr='Sep',dt_days=5,
                    also_return_relative_dts=False,
                    relative_dts_inplace=True,
                    verbose=False):

    if not hasattr(dt_days,'__len__'):
        dt = timedelta(days=dt_days)
        dt_lower = dt
        dt_upper = dt
        dtdaystr = "+/- {:.2f} days".format(dt_days)
    else:
        assert len(dt_days) == 2
        dt_lower = timedelta(dats=dt_days[0])
        dt_upper = timedelta(dats=dt_days[1])
        dtdaystr = "+ {:.2f} d/- {:.2f} days".format(dt_days[1],dt_days[0])

    if verbose:
        print("get_season_inds: df has {:8d} inds, requested season is {:5s} {:s}".
              format(df.shape[0],seasonStr,dtdaystr))


    # Initialize all indices to FALSE
    indices = df.index < df.index.min()

    if also_return_relative_dts:
        if relative_dts_inplace:
            if 'reldoy' not in list(df.columns):
                df['reldoy'] = -99999

        else:
            relDTs = df.index - df.index[0]
            print("Does it work to fill the whole column up with NaTs?")
            breakpoint()
            relDTs[:] = pd.NaT

    seasonDict = load_fancy_season_dict(y0=df.index.year.min()-1,
                                        y1=df.index.year.max()+1,
                                        return_doyDict=False,
                                        drop_timezone_info=True)

    for yr,seasonDT in seasonDict[seasonStr].items():
        # print(seasonDT)
    
        these = (df.index >= (seasonDT-dt_lower)) & (df.index <= (seasonDT+dt_upper))
        
        nHere = np.where(these)[0].size
        if nHere > 0:
            # print("Dodat med {:d}".format(nHere))
            if verbose:
                print("Got {:5d} inds for {:4d} {:5s}".format(nHere,yr,seasonStr))

            indices[these] = True

            if also_return_relative_dts:
                if relative_dts_inplace:
                    
                    df.loc[these,'reldoy'] = (df[these].index - seasonDT).total_seconds().values/(24*3600.)

                else:
                    relDTs[these] = df.index[these] - seasonDT

    if also_return_relative_dts:
        if relative_dts_inplace:
            return indices
        else:
            return indices,relDTs
    else:
        return indices


def load_fancy_season_dict(y0=1900,
                           y1=2051,
                           return_doyDict=False,
                           drop_timezone_info=True):
    """
    seasonDict = load_fancy_season_dict(y0=1990,y1=2020)
    Get a dictionary of March and September equinox times + June and December solstice times

    Based on: https://rhodesmill.org/skyfield/almanac.html
    """

    ########################################
    # COMPARE WITH load_season_dict:
    ########################################
    # seasonDict = load_season_dict(load_extended=bigMode)
    # seasonDict2 = load_fancy_season_dict(y0=2010,
    #                                      y1=2020,
    #                                      return_doyDict=False,
    #                                      drop_timezone_info=True)
    #
    # checkseason = 'Sep'
    # for yr in range(2010,2020+1):
    #     print(yr,seasonDict2[checkseason][yr],np.abs(seasonDict2[checkseason][yr]-seasonDict[checkseason][str(yr)]))

    doOldWay = not drop_timezone_info

    if return_doyDict:
        typeStr = 'doy'
    else:
        typeStr = 'datetime'

    print("Getting {:s}Dict for {:d}-{:d} ...".format(typeStr, y0, y1))


    seasonDict = load_season_dict_local(y0=y0,
                                        y1=y1,
                                        return_doyDict=return_doyDict)

    if not doOldWay:
        return seasonDict

    # OLDSTUFF

    # OLD SKYFIELD LOADER
    # ts = api.load.timescale()
    # e = api.load('de421.bsp')

    # NEW SKYFIELD LOADER

    from skyfield import api
    from skyfield import almanac

    if return_doyDict:
        from hatch_python_utils.time_tools import datetime_to_doy

    from hatch_python_utils.hatch_utils import get_basedir
    basedir, isColtrane = get_basedir(verbose=True)
    
    skyfielddir = basedir + 'database/skyfield/'
    load = api.Loader(skyfielddir)
    ts = load.timescale()
    e = load('de421.bsp')

    # Note that almanac computation can be slow and expensive. To determine the moment of sunrise, for example, Skyfield has to search back and forth through time asking for the altitude of the Sun over and over until it finally works out the moment at which it crests the horizon.

    # Create a start time and an end time to ask for all of the equinoxes and solstices that fall in between.
    years = [year for year in range(y0, y1+1)]
    seasonDict = dict(Mar={}, Jun={}, Sep={}, Dec={})

    if return_doyDict:
        typeStr = 'doy'
    else:
        typeStr = 'datetime'

    print("Getting {:s}Dict for {:d}-{:d} ...".format(typeStr, y0, y1))

    # for year in years:

    #     # print(year)

    #     t0 = ts.utc(year, 1, 1)
    #     t1 = ts.utc(year, 12, 31)
    #     t, y = almanac.find_discrete(t0, t1, almanac.seasons(e))

    #     # for yi, ti in zip(y, t):
    #     #     print(yi, almanac.SEASON_EVENTS[yi], ti.utc_iso(' '))

    #     if return_doyDict:
    #         doys = datetime_to_doy([t[0].utc_datetime(),
    #                                 t[1].utc_datetime(),
    #                                 t[2].utc_datetime(),
    #                                 t[3].utc_datetime()])

    #         seasonDict['Mar'][year] = doys[0]
    #         seasonDict['Jun'][year] = doys[1]
    #         seasonDict['Sep'][year] = doys[2]
    #         seasonDict['Dec'][year] = doys[3]

    #     else:
    #         if drop_timezone_info:
    #             seasonDict['Mar'][year] = t[0].utc_datetime().replace(
    #                 tzinfo=None)
    #             seasonDict['Jun'][year] = t[1].utc_datetime().replace(
    #                 tzinfo=None)
    #             seasonDict['Sep'][year] = t[2].utc_datetime().replace(
    #                 tzinfo=None)
    #             seasonDict['Dec'][year] = t[3].utc_datetime().replace(
    #                 tzinfo=None)
    #         else:
    #             seasonDict['Mar'][year] = t[0].utc_datetime()
    #             seasonDict['Jun'][year] = t[1].utc_datetime()
    #             seasonDict['Sep'][year] = t[2].utc_datetime()
    #             seasonDict['Dec'][year] = t[3].utc_datetime()

    return seasonDict


def load_season_dict_local(y0=1900,
                           y1=2051,
                           return_doyDict=False):

    assert (y0 >= 1900) & (y1 <= 2051)

    if return_doyDict:
        from hatch_python_utils.time_tools import datetime_to_doy

    TEST_FILENAME = os.path.join(os.path.dirname(__file__), '../dataz/earth_seasons_1900-2051.csv')

    file1 = open(TEST_FILENAME, 'r') 
    Lines = file1.readlines() 

    strftimefmt = "%Y-%m-%d %H:%M:%S.%f"
    seasonDict2 = dict(Mar=dict(),Jun=dict(),Sep=dict(),Dec=dict())
    for line in Lines:

        yr = int(line[0:4])

        if yr < y0:
            continue

        if yr > y1:
            continue

        seps = line.split(",")
        # print(yr)

        seasonDict2['Mar'][yr] = datetime.strptime(seps[0],strftimefmt)
        seasonDict2['Jun'][yr] = datetime.strptime(seps[1],strftimefmt)
        seasonDict2['Sep'][yr] = datetime.strptime(seps[2],strftimefmt)
        seasonDict2['Dec'][yr] = datetime.strptime(seps[3].replace('\n',''),strftimefmt)
    

        if return_doyDict:
            doys = datetime_to_doy([seasonDict2['Mar'][yr],
                                    seasonDict2['Jun'][yr],
                                    seasonDict2['Sep'][yr],
                                    seasonDict2['Dec'][yr]])

            seasonDict2['Mar'][yr] = doys[0]
            seasonDict2['Jun'][yr] = doys[1]
            seasonDict2['Sep'][yr] = doys[2]
            seasonDict2['Dec'][yr] = doys[3]

    return seasonDict2


def load_season_dict(load_extended=False):

    # tides
    if load_extended:
        tides = [['20010320 13:31', '20010621 07:38',
                  '20010922 23:05', '20011221 19:22'],
                 ['20020320 19:16', '20020621 13:25',
                  '20020923 04:56', '20021222 01:15'],
                 ['20030321 01:00', '20030621 19:11',
                  '20030923 10:47', '20031222 07:04'],
                 ['20040320 06:49', '20040621 00:57',
                  '20040922 16:30', '20041221 12:42'],
                 ['20050320 12:34', '20050621 06:46',
                  '20050922 22:23', '20051221 18:35'],
                 ['20060320 18:25', '20060621 12:26',
                  '20060923 04:04', '20061222 00:22'],
                 ['20070321 00:07', '20070621 18:06',
                  '20070923 09:51', '20071222 06:08'],
                 ['20080320 05:49', '20080621 00:00',
                  '20080922 15:45', '20081221 12:04'],
                 ['20090320 11:44', '20090621 05:45',
                  '20090922 21:18', '20091221 17:47'],
                 ['20100320 17:32', '20100621 11:28',
                  '20100923 03:09', '20101221 23:38'],
                 ['20100320 17:32', '20100621 11:28',
                  '20100923 03:09', '20101221 23:38'],
                 ['20110320 23:21', '20110621 17:16',
                  '20110923 09:04', '20111222 05:30'],
                 ['20120320 05:14', '20120620 23:09',
                  '20120922 14:49', '20121221 11:12'],
                 ['20130320 11:02', '20130621 05:04',
                  '20130922 20:44', '20131221 17:11'],
                 ['20140320 16:57', '20140621 10:51',
                  '20140923 02:29', '20141221 23:03'],
                 ['20150320 22:45', '20150621 16:38',
                  '20150923 08:20', '20151222 04:48'],
                 ['20160320 04:30', '20160620 22:34',
                  '20160922 14:21', '20161221 10:44'],
                 ['20170320 10:28', '20170621 04:24',
                  '20170922 20:02', '20171221 16:28'],
                 ['20180320 16:15', '20180621 10:07',
                  '20180923 01:54', '20181221 22:23'],
                 ['20190320 21:58', '20190621 15:54',
                  '20190923 07:50', '20191222 04:19'],
                 ['20200320 03:50', '20200620 21:44',
                  '20200922 13:31', '20201221 10:02']]
    else:
        tides = [['20100320 17:32', '20100621 11:28',
                  '20100923 03:09', '20101221 23:38'],
                 ['20110320 23:21', '20110621 17:16',
                  '20110923 09:04', '20111222 05:30'],
                 ['20120320 05:14', '20120620 23:09',
                  '20120922 14:49', '20121221 11:12'],
                 ['20130320 11:02', '20130621 05:04',
                  '20130922 20:44', '20131221 17:11'],
                 ['20140320 16:57', '20140621 10:51',
                  '20140923 02:29', '20141221 23:03'],
                 ['20150320 22:45', '20150621 16:38',
                  '20150923 08:20', '20151222 04:48'],
                 ['20160320 04:30', '20160620 22:34',
                  '20160922 14:21', '20161221 10:44'],
                 ['20170320 10:28', '20170621 04:24',
                  '20170922 20:02', '20171221 16:28'],
                 ['20180320 16:15', '20180621 10:07',
                  '20180923 01:54', '20181221 22:23'],
                 ['20190320 21:58', '20190621 15:54',
                  '20190923 07:50', '20191222 04:19'],
                 ['20200320 03:50', '20200620 21:44',
                  '20200922 13:31', '20201221 10:02']]

    seasonDict = dict(Mar={}, Jun={}, Sep={}, Dec={})
    for i, key in enumerate(seasonDict.keys()):
        for tide in tides:
            tmpdate = datetime.strptime(tide[i], "%Y%m%d %H:%M")
            seasonDict[key][str(tmpdate.year)] = tmpdate

    print("Adding 'spring', 'summer', 'fall', and 'winter' keys for backwards-compatibility")
    seasonDict['spring'] = seasonDict['Mar']
    seasonDict['summer'] = seasonDict['Jun']
    seasonDict['fall'] = seasonDict['Sep']
    seasonDict['winter'] = seasonDict['Dec']

    return seasonDict


def get_scaled_season_parameter(timestamps,
                                verbose=False):
    """
    tau = get_scaled_season_parameter(timestamps)

    Scales timestamps such that 
    tau = 0: March equinox
    tau = 1: June solstice
    tau = 2: September equinox
    tau = 3: December solstice

    Parameter
    =========

    timestamps : pandas.DatetimeIndex (or perhaps numpy.datetime64 array?)
    """
    
    # Need pd.DatetimeIndex
    if not isinstance(timestamps,pd.DatetimeIndex):
        timestamps = pd.DatetimeIndex(timestamps)

    reversesesongdict = dict(Mar='Jun',
                             Jun='Sep',
                             Sep='Dec',
                             Dec='Mar')
    
    y0 = timestamps.year.min()-1
    y1 = timestamps.year.max()+1

    seasonDict = load_fancy_season_dict(y0=y0,
                                        y1=y1,
                                        return_doyDict=False,
                                        drop_timezone_info=True)

    # Scale season thing
    sesongting = dict(zip(list(seasonDict.keys()),np.arange(0,4)))
    
    tau = np.zeros(timestamps.size,dtype=np.float64)*np.nan
    
    checkyears = list(seasonDict['Mar'].keys())
    
    count = 0
    for chkyr in checkyears:
        for sesong,offset in sesongting.items():
    
            sesongtid = seasonDict[sesong][chkyr]
            if sesong == 'Dec':
                nestesesongtid = seasonDict['Mar'][np.clip(chkyr+1,np.min(checkyears),np.max(checkyears))]
            else:
                nestesesongtid = seasonDict[reversesesongdict[sesong]][chkyr]
    
            dtsesong = nestesesongtid-sesongtid
            if verbose:
                print(count,chkyr,sesong,offset,
                      sesongtid.strftime("%Y-%m-%d %H:%M:%S"),
                      nestesesongtid.strftime("%Y-%m-%d %H:%M:%S"),
                      dtsesong)
    
            haveinds = (timestamps >= sesongtid) & (timestamps < nestesesongtid)
            nHaveInds = haveinds.sum()
    
            if nHaveInds == 0:
                if verbose:
                    print("None here!")
                continue
    
            if verbose:
                print("Got {:d} inds".format(nHaveInds))
    
            tau[haveinds] = (timestamps[haveinds]-sesongtid).total_seconds().values/dtsesong.total_seconds()+offset

            count += 1

        #     if count >= 10:
        #         break
    
        # if count >= 10:
        #     break
        #     break
        # break

    return tau


def add_scaled_season_column(df,
                             data_latitudecol='alat110',
                             verbose=False):

    df['scaledtidoffset'] = 0.

    df['scaledtidoffset'] = get_scaled_season_parameter(df.index,
                                                        verbose=verbose)

    df['localscaledtidoffset'] = df['scaledtidoffset']
    if data_latitudecol in df.columns:
        shinds = df[data_latitudecol] < 0
        df.loc[shinds,'localscaledtidoffset'] = (df.loc[shinds,'localscaledtidoffset'] + 2) % 4
    else:
        print(f"Column '{data_latitudecol}' does not exist in provided df! Not calculating 'localscaledtidoffset'")
