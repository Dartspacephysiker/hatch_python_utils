import pandas as pd
import numpy as np
from datetime import datetime,timedelta


def earthSunDist(doy):
    """
    https://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year

    This is an approximate expression. Term by term,

    1       : The mean distance between the Earth and the Sun is about one
    astronomical unit.

    0.01672 : This is the eccentricity of the Earth's orbit about the Sun.

    cos     : This is of course the cosine function, but with argument in 
              degrees rather than radians.

    0.9856  : This is 360/365.256363, where 360 is the number of degrees in a
              full rotation and 365.256363 is the length of a sidereal year, in
              mean solar days.

    doy     : This is the day number of the year. Since this is an approximate
              expression, whether one starts with the first of January being
              zero or one is irrelevant.

    4       : The Earth currently reaches perihelion between the fourth and sixth of
              January, depending on the year.

    RETURNS : Earth-Sun distance as a fraction of 1 AU (= 149,597,871 km)

    *********************
    Example
    *********************
    import numpy as np
    from hatch_python_utils import earth as he
    figD = plt.figure()
    axD = plt.subplot(1,1,1)

    doys = np.arange(1,365)
    junk = axD.plot(doys,he.earthSunDist(doys))
    """
    return 1-0.01672*np.cos(np.deg2rad(0.9856*(doy-4)))


def sphDist(lat1, mltOrLon1, lat2, mltOrLon2,
            mltMode=True):
    """
    Great-circle distance
    lat1: Latitude of 1st point in degrees
    lat2: Latitude of 2nd point in degrees
    mltOrLon1: Magnetic local time of 1st point in hours (so between 0 and 24)
    mltOrLon2: Magnetic local time of 2nd point in hours
    """
    lat1R = np.deg2rad(lat1)
    lat2R = np.deg2rad(lat2)
    if mltMode:
        lon1R = np.deg2rad(mltOrLon1*15.)
        lon2R = np.deg2rad(mltOrLon2*15.)
    else:
        lon1R = np.deg2rad(mltOrLon1)
        lon2R = np.deg2rad(mltOrLon2)

    return np.rad2deg(np.arccos(np.sin(lat1R)*np.sin(lat2R)+np.cos(lat1R)*np.cos(lat2R)*np.cos(np.abs(lon1R-lon2R))))

# def sphDist(lat1, mlt1, lat2, mlt2):
#     """
#     Great-circle distance
#     lat1: Latitude of 1st point in degrees
#     lat2: Latitude of 2nd point in degrees
#     mlt1: Magnetic local time of 1st point in hours (so between 0 and 24)
#     mlt2: Magnetic local time of 2nd point in hours
#     """
#     lat1R = np.deg2rad(lat1)
#     mlt1R = np.deg2rad(mlt1*15.)
#     lat2R = np.deg2rad(lat2)
#     mlt2R = np.deg2rad(mlt2*15.)

#     return np.rad2deg(np.arccos(np.sin(lat1R)*np.sin(lat2R)+np.cos(lat1R)*np.cos(lat2R)*np.cos(np.abs(mlt1R-mlt2R))))

# # geodeticheight2geocentricR(lat, height)


#     gdlatJ, altJ, XIGRF, ZIGRF = geodesy.geoc2geod(90.-geolat,
#                                                    dfMInterp["Radius"].values/1000.,
#                                                    -igrf.north.values, -igrf.down.values)
#     igrfMag = np.sqrt(igrf.east.values**2.+igrf.north.values **
#                       2. + igrf.down.values**2.)

#     return

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


    from skyfield import api
    from skyfield import almanac

    if return_doyDict:
        from hatch_python_utils.time_tools import datetime_to_doy

    ts = api.load.timescale()
    e = api.load('de421.bsp')

    # Note that almanac computation can be slow and expensive. To determine the moment of sunrise, for example, Skyfield has to search back and forth through time asking for the altitude of the Sun over and over until it finally works out the moment at which it crests the horizon.

    # Create a start time and an end time to ask for all of the equinoxes and solstices that fall in between.
    years = [year for year in range(y0, y1+1)]
    seasonDict = dict(Mar={}, Jun={}, Sep={}, Dec={})

    if return_doyDict:
        typeStr = 'doy'
    else:
        typeStr = 'datetime'

    print("Getting {:s}Dict for {:d}-{:d} ...".format(typeStr, y0, y1))

    for year in years:

        # print(year)

        t0 = ts.utc(year, 1, 1)
        t1 = ts.utc(year, 12, 31)
        t, y = almanac.find_discrete(t0, t1, almanac.seasons(e))

        # for yi, ti in zip(y, t):
        #     print(yi, almanac.SEASON_EVENTS[yi], ti.utc_iso(' '))

        if return_doyDict:
            doys = datetime_to_doy([t[0].utc_datetime(),
                                    t[1].utc_datetime(),
                                    t[2].utc_datetime(),
                                    t[3].utc_datetime()])

            seasonDict['Mar'][year] = doys[0]
            seasonDict['Jun'][year] = doys[1]
            seasonDict['Sep'][year] = doys[2]
            seasonDict['Dec'][year] = doys[3]

        else:
            if drop_timezone_info:
                seasonDict['Mar'][year] = t[0].utc_datetime().replace(
                    tzinfo=None)
                seasonDict['Jun'][year] = t[1].utc_datetime().replace(
                    tzinfo=None)
                seasonDict['Sep'][year] = t[2].utc_datetime().replace(
                    tzinfo=None)
                seasonDict['Dec'][year] = t[3].utc_datetime().replace(
                    tzinfo=None)
            else:
                seasonDict['Mar'][year] = t[0].utc_datetime()
                seasonDict['Jun'][year] = t[1].utc_datetime()
                seasonDict['Sep'][year] = t[2].utc_datetime()
                seasonDict['Dec'][year] = t[3].utc_datetime()

    return seasonDict


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


def cheap_LT_calc(dts,gclons,
                  return_dts_too=False,
                  verbose=False):

    # print("Convert input dts to pandas datetime index!")
    if not hasattr(dts,'__iter__'):
        dts = pd.DatetimeIndex([dts])

    elif not hasattr(dts,'hour'):
        dts = pd.DatetimeIndex(dts)

    if not hasattr(gclons,'__iter__'):
        gclons = np.array(gclons)

    gclons = (gclons+360) % 360

    if verbose:
        if gclons.size ==1:
            # relstr = "ahead of"
            reltogreenwich = gclons/15.

            relstr = "ahead of" if (gclons <= 180) else "behind"
            if reltogreenwich > 12:
                reltogreenwich -= 24

            print("Longitude {:.2f} is {:.2f} hours {:s} Greenwich!".format(gclons,reltogreenwich,relstr))

    midnightlongitude = -15*(dts.hour.values+dts.minute.values/60+dts.second.values/3600.)
    midnightlongitude = (midnightlongitude + 360) % 360

    LTs = (((gclons-midnightlongitude) + 360) % 360)/15
    if return_dts_too:
        return LTs, dts, midnightlongitude, gclons
    else:
        return LTs

def get_noon_longitude(dts,verbose=False):
    """
    """
    # A test:
    # import datetime
    # import pandas as pd
    # from pytt.earth.sunlight import sza
    # 
    # marequinox = datetime.datetime(2015, 3, 20, 22, 45, 9, 340000)
    # junsolstice = datetime.datetime(2015, 6, 21, 16, 37, 55, 813000)
    # refdate = junsolstice
    # refdate = marequinox
    # dts = pd.date_range(start=datetime(refdate.year,refdate.month,refdate.day,0),
    #                     end=datetime(refdate.year,refdate.month,refdate.day,23),
    #                     freq='1h')
    # # FOLLOWING SHOULD ALL BE AROUND 23.44 IF JUNSOLSTICE, 0 IF MAREQUINOX
    # print(sza(np.zeros(dts.size),get_noon_longitude(dts),dts))

    if not hasattr(dts,'__iter__'):
        dts = pd.DatetimeIndex([dts])

    elif not hasattr(dts,'hour'):
        dts = pd.DatetimeIndex(dts)

    fracHour = dts.hour.values+dts.minute.values/60+dts.second.values/3600.

    assert not any((fracHour < 0) | (fracHour > 24))

    fracHour[fracHour > 12] -= 24

    fracHour *= 15

    if verbose:
        print("Min fracHour: {:.2f}".format(np.min(fracHour)))
        print("Max fracHour: {:.2f}".format(np.max(fracHour)))

    return 180 - fracHour 

