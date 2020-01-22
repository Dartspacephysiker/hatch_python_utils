import pandas as pd
import numpy as np
from datetime import datetime,timedelta
from .geodesy import get_max_sza
from pytt.earth.sunlight import sza


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
    """
    # Basert på dette: https://rhodesmill.org/skyfield/almanac.html

    from skyfield import api
    from skyfield import almanac

    if return_doyDict:
        from hatch_python_utils.time_tools import datetime_to_doy

    ts = api.load.timescale()
    e = api.load('de421.bsp')

    # Note that almanac computation can be slow and expensive. To determine the moment of sunrise, for example, Skyfield has to search back and forth through time asking for the altitude of the Sun over and over until it finally works out the moment at which it crests the horizon.

    # Create a start time and an end time to ask for all of the equinoxes and solstices that fall in between.
    years = [year for year in range(y0, y1)]
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

    seasonDict = dict(spring={}, summer={}, fall={}, winter={})
    for i, key in enumerate(seasonDict.keys()):
        for tide in tides:
            tmpdate = datetime.strptime(tide[i], "%Y%m%d %H:%M")
            seasonDict[key][str(tmpdate.year)] = tmpdate

    return seasonDict

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

def get_t_in_darkness(alt,glat,glon,dates,
                      tol_deg=0.15,
                      verbose=True,
                      dodiagnosticprint=False):
    """
    For a given timestamp, altitude and GEOGRAPHIC (or is it geodetic??? Ask Kalle what sza uses!) lat and lon, calculate how many seconds
    this point has been in darkness assuming no refraction of light and a perfectly spherical earth.
    
    alt           : Geodetic altitude       (km)
    glat          : Geographic(?) latitude  (degrees)
    glon          : Geographic(?) longitude (degrees)
    dates         : Timestamps              (datetime, pandas DatetimeIndex, etc.)
    tol_deg       : Fudge factor            (degrees).
                    I find that a fraction of a degree (specifically 0.15) seems to do the job
    """

    #So how do we do this?
    #
    ##########
    # 1. Get solar zenith angle (sza)
    #
    ##########
    # 2. Get maximum sza at which sun is visible for given altitude (assumes spherical Earth!)
    #    During this step we find out which points lie in darkness (taking stock of 'tol_deg' fudge factor).
    #    We mark points that are already sunlit as such, and give them t_in_darkness = 0.
    #
    ##########
    # 3. Get SZA for each altitude, latitude, and NOON longitude; see if each latitude is CURRENTLY sunlit
    #    at the local-noon longitude.
    #
    ##########
    # 4. For those alt/tstamp/lat point for which the local-noon longitude is NOT sunlit, iteratively shift
    #    tHadSun back by one day until we find that there is sunshine at local-noon longitude.
    #
    #    After this step, tHadSun will be an array of timestamps for which the sun is visible at the given 
    #    alt/lat and calculated NOON-longitude
    #
    # 4a. A check – all noonszas must be less than their corresponding maximum sza (i.e., all noon-szas must 
    #     now be sunlit) before step 5.
    #
    ##########
    # 5. Calculate the time shift (which, as a result of step 4, is at most 23.9999 hours) needed to put each 
    #    alt/lat/lon pair at noon. Subtract this time shift from tHadSun so that tHadSun corresponds to the
    #    last day on which this alt/lat/lon pair was in sunlight at local noon.
    #
    # 5a. Do some fudging here.  This fudging is necessary because, for a given latitude, the minimum sza
    # obtained over the course of a day changes.
    #
    # FREE EXAMPLE TO EXPLAIN WHY WE HAVE TO FUDGE: 
    # Consider an alt/lat/lon that is in darkness at, say, 03:00 local time. Now imagine following this point
    # along a line of constant altitude and latitude for the given timestamp (i.e. vary longitude ONLY and
    # hold all else constant) until we reach the longitude corresponding to local noon. Suppose that this
    # alt/lat is sunlit at local noon. Now, because the point we are actually interested currently lies at
    # 03:00 local time, we have to subtract 15 hours from the current time stamp to put our alt/lat/lon at
    # local noon. But this particular alt/lat pair may not have been sunlit at local noon 15 hours before the
    # present time! The necessary fudge factor appears to be of order 0.1 degrees.
    #
    ##########
    # 6. After shifting timestamps to put this longitude at local noon, all alt/lat/lon pairs are now sunlit
    # (up to the fudge factor 'tol_deg'). Now we just increment each timestamp until the alt/lat/lon pair
    # falls in darkness. The stepping in time goes by hours, then minutes, then seconds. 
    
    # 1. Get sza
    origsza = sza(glat,glon,dates)

    # 2. Get thresh SZA for given altitude
    maxsza = get_max_sza(alt)
    
    alreadysunlit = origsza <= (maxsza+np.abs(tol_deg))
    origDark = origsza > maxsza
    nOrigDarkPts = np.sum(origDark)

    if verbose:
        print("{:d} location/time pairs in darkness ({:d} in light)!".
              format(np.sum(~alreadysunlit),
                     np.sum(alreadysunlit)))

    if nOrigDarkPts == 0:
        return np.zeros(glat.size)

    # 3. Get SZA for each altitude, latitude, and NOON longitude to see if each point is/may have been in sunlight
    noonsza = sza(glat,
                  get_noon_longitude(dates),
                  dates)
    
    tHadSun = dates.copy()
    
    stillDark = noonsza > maxsza
    fixed = noonsza <= maxsza       # Keep track of which points need updating
    nStillDarkPts = stillDark.sum()
    
    if verbose:
        print("{:d} of these latitudes are in darkness at local-noon longitude ({:d} in light)!".
              format(nStillDarkPts,fixed.sum()))
    
    # 4. For each point, shift tHadSun back by one day until we find that there is sunshine at local-noon longitude.
    # After this step, tHadSun will be an array of timestamps for which the sun is visible for the given latitude
    # and altitude, and calculated NOON-longitude 
    
    daysback = 1
    totalNoonAdjusted = 0
    while nStillDarkPts > 0:
    
        thistdelta = timedelta(days=daysback)
    
        # DIAG
        if dodiagnosticprint:
            print(tHadSun[stillDark]-thistdelta)
    
        noonsza[stillDark] = sza(glat[stillDark],
                                 get_noon_longitude(tHadSun[stillDark]),
                                 tHadSun[stillDark]- thistdelta)
    
        # Calculate who is still in darkness, N stilldark points
        stillDark = noonsza > maxsza
        fixme = ~(stillDark) & ~(fixed) & ~(alreadysunlit)
        nFixme = np.sum(fixme)
    
        if nFixme > 0:
            
            totalNoonAdjusted += nFixme
            if dodiagnosticprint:
                print("Moving {:d} timestamps !".format(nFixme))
    
            tHadSun[fixme] = tHadSun[fixme]-thistdelta
    
            fixed[fixme] = True
    
        nStillDarkPts = stillDark.sum()
        daysback += 1
    
    if verbose:
        print("Moved {:d} timestamps back".format(totalNoonAdjusted))

    # 4a. A check – all noonszas should be less than their corresponding maxsza
    noonsza = sza(glat,get_noon_longitude(tHadSun),tHadSun)
    assert all(noonsza[~alreadysunlit] <= maxsza[~alreadysunlit])
    
    # 5. Calculate the time shift (which, as a result of step 4, is at most 23.9999 hours) needed to put each alt/lat/lon pair at noon.
    #    Subtract this time shift from tHadSun so that tHadSun corresponds to the last day on which this alt/lat/lon pair was in sunlight at local noon.
    shiftHours = cheap_LT_calc(tHadSun,
                               glon,
                               return_dts_too=False,
                               verbose=True)

    shiftHours = shiftHours - 12
    shiftHours[shiftHours < 0] += 24

    shiftHours[alreadysunlit] = 0

    timedeltas = pd.TimedeltaIndex(
        data=shiftHours*3600,
        unit='s').to_pytimedelta()
    
    # Don't believe this is noon for these guys? Try it:
    # print((cheap_LT_calc(tHadSun-timedeltas,
    #                      glon,
    #                      return_dts_too=False,
    #                      verbose=True)).describe())

    # RESET STILLDARK TO ACTUALDARKS
    testsza = origsza.copy()
    stillDark = origDark.copy()
    nStillDarkPts = np.sum(origDark)
    
    testsza[stillDark] = sza(glat[stillDark],
                             glon[stillDark],
                             tHadSun[stillDark]-timedeltas[stillDark])
    
    # assert all(testsza[stillDark] < maxsza[stillDark])

    # 5a. Do some fudging here.  This fudging is necessary because, for a given latitude, the minimum sza obtained over the course of a
    # day changes.
    # 
    if not all(testsza[stillDark] <= maxsza[stillDark]):

        diff = (testsza[stillDark] - maxsza[stillDark])

        if diff.max() > tol_deg:
            print("Warning! error of more than {:.2f} deg in darknesscalc!".format(tol_deg))
            breakpoint()

        badHeads = np.where(diff > 0)[0]
        badHeads = np.where(stillDark)[0][badHeads]

        maxsza[badHeads] += diff[diff > 0]
        #     print("N badheads: {:d}. Rotate 'em back a day".format(badHeads.size))
        #     tHadSun[badHeads] -= timedelta(days=1)

    dodat = stillDark & ~alreadysunlit
    tHadSun[dodat] = tHadSun[dodat]-timedeltas[dodat]
    
    # 6. After shifting timestamps to put this longitude at local noon, all alt/lat/lon pairs are now sunlit
    # (up to the fudge factor 'tol_deg'). Now we just increment each timestamp until the alt/lat/lon pair falls in darkness.
    # The stepping in time goes by hours, then minutes, then seconds. 

    #7. If some need to be fixed, go back each minute until we're where we need to be 
    # No er det berre å trylle tiden framover for å finne tidspunktet hvor mørket slår
    # tHadSun[stillDark] is filled with points that are now LIGHT
    
    daystep = 0
    hourstep = 1
    minutestep = 0
    secondstep = 0
    stepcount = 0
    
    omgangtype = 'hours'
    thistdelta = timedelta(days=daystep,hours=hourstep,minutes=minutestep,seconds=secondstep)
    
    haveSteppedDays = True
    haveSteppedHours = False
    haveSteppedMinutes = False
    haveSteppedSeconds = False
    
    testsza = origsza.copy()
    nStillLightPts = np.sum(origDark)
    stillLight = origDark.copy()
    
    while (nStillLightPts > 0) and (not haveSteppedSeconds):
    
        oldtestsza = testsza.copy()

        # Get sza for darkpoints given this timedelta
        testsza[stillLight] = sza(glat[stillLight],
                                 glon[stillLight],
                                 tHadSun[stillLight]+thistdelta)
    
        # Calculate who is still in light, N stillLight points
        # stillLight = (testsza < (maxsza-tol_deg)) & ~alreadysunlit
        stillLight = (testsza <= maxsza) & ~alreadysunlit
        nStillLightPts = stillLight.sum()

        if (omgangtype == 'hours') & (stepcount >= 24):
            print("Bogusness")
            breakpoint()

        if stepcount > 0:
            if np.where(testsza[stillLight] < oldtestsza[stillLight])[0].size > 0:
                print("BAD!")

        if dodiagnosticprint:
            print("Adjusting tstamp for {:d} points!".format(np.sum(stillLight)))
    
        # Update timestamps for those that are still light, even with this time adjustment
        tHadSun[stillLight] += thistdelta
    
        if np.where(tHadSun > dates)[0].size > 0:
            print("Bogus!")
            breakpoint()

        print("nStillLight: ",nStillLightPts)

        # if nStillLightPts == 1:
        #     breakpoint()

        if nStillLightPts == 0:
            if dodiagnosticprint:
                print("No more lights for {:s} omgang!".format(omgangtype))
    
            if not haveSteppedDays:
                haveSteppedDays = True
    
                daystep = 0
                hourstep = 1
                minutestep = 0
                secondstep = 0
                stepcount = 0
                thistdelta = timedelta(days=daystep,
                                       hours=hourstep,
                                       minutes=minutestep,
                                       seconds=secondstep)
                omgangtype = 'hours'
    
                testsza = origsza.copy()
                nStillLightPts = np.sum(origDark)
                stillLight = origDark.copy()
    
            elif not haveSteppedHours:
                haveSteppedHours = True
    
                daystep = 0
                hourstep = 0
                minutestep = 1
                secondstep = 0
                stepcount = 0
                thistdelta = timedelta(days=daystep,
                                       hours=hourstep,
                                       minutes=minutestep,
                                       seconds=secondstep)
                omgangtype = 'minutes'
    
                testsza = origsza.copy()
                nStillLightPts = np.sum(origDark)
                stillLight = origDark.copy()
    
            elif not haveSteppedMinutes:
                haveSteppedMinutes = True
    
                daystep = 0
                hourstep = 0
                minutestep = 0
                secondstep = 1
                stepcount = 0
                thistdelta = timedelta(days=daystep,
                                       hours=hourstep,
                                       minutes=minutestep,
                                       seconds=secondstep)
                omgangtype = 'seconds'
    
                testsza = origsza.copy()
                nStillLightPts = np.sum(origDark)
                stillLight = origDark.copy()
    
            elif not haveSteppedSeconds:
                haveSteppedSeconds = True
    
        stepcount += 1
        # if dodiagnosticprint:
        print("{:4d} {:s} steps".format(stepcount,omgangtype))
    
    # sza(glat,
    #     glon,
    #     tHadSun)

    # finalcheck
    if not all(sza(glat,
                   glon,
                   tHadSun) <= maxsza):
        breakpoint()

    return pd.TimedeltaIndex(dates-tHadSun).total_seconds()
