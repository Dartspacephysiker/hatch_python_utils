#2020/02/28
# def sunlight():

import pandas as pd
import numpy as np
from datetime import datetime,timedelta

from pytt.earth.sunlight import sza


def get_max_sza(z,
                R=6371.):
    """
    z is altitude in km
    # R = 6371. # Earth radius (default)
    """

    # assert not hasattr(R,'__iter__')
    zIsArray = hasattr(z,'size')

    if not zIsArray:
        z = np.array([z])

    h = z + R
    max_sza = 180.-np.rad2deg(np.arctan2(np.sqrt(1-(R/h)**2.), (h/R-R/h)))

    # if zIsArray:
        
    max_sza[np.isclose(z,0)] = 90.
        
    max_sza[z < 0] = np.nan
        
    # else:
        
    #     if np.isclose(z,0):
    #         max_sza = 90.
    #     elif (z < 0):
    #         max_sza = np.nan

    if zIsArray:
        return max_sza
    else:
        return max_sza[0]


def get_daily_sza_stats(gdlat,gdlon,date,freq='15min',
                        degrees=True,
                        statfunctiondict=None,
                        return_dataframe=True,
                        add_sunlit_darkness_stats=True,
                        sunlit_darkness__alt=500):

    assert (date.hour + date.minute/60. + date.second/3600.) == 0

    if statfunctiondict is not None:

        stats = statfunctiondict

    else:

        stats = dict(mean=np.mean,
                 median=np.median,
                 std=np.std,
                 min=np.min,
                 max=np.max,
                 Q1=lambda x: np.quantile(x,0.25),
                 Q3=lambda x: np.quantile(x,0.75))

    Nlatlon = gdlat.size

    if add_sunlit_darkness_stats:
        max_sza = get_max_sza(sunlit_darkness__alt,
                              R=6371.)

    tmpdates = pd.date_range(start=date,
                             end=date+timedelta(days=1),
                             freq=freq,
                             closed='left')

    szas = np.zeros((tmpdates.size,Nlatlon))

    for idate,tmpdate in enumerate(tmpdates):
        szas[idate,:] = sza(gdlat, gdlon, tmpdate, 
                            degrees=degrees)

    szastats = dict()    

    for key,func in stats.items():
        # print(key,func)
        szastats[key] = func(szas)

    if add_sunlit_darkness_stats:
        
        sunlit = szas < max_sza

        sunlitfracs = np.sum(sunlit,axis=1)/Nlatlon
        szastats['sunlit_maxfrac'] = np.max(sunlitfracs)
        szastats['sunlit_minfrac'] = np.min(sunlitfracs)
        szastats['sunlit_meanfrac'] = np.mean(sunlitfracs)

    if return_dataframe:
        return pd.DataFrame(szastats,index=pd.DatetimeIndex([date]))
    else:
        return szastats	


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
