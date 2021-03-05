import pandas as pd
import numpy as np
from datetime import datetime,timedelta
from pysymmetry.sunlight import subsol

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

def get_localtime(dts, glons,
                  return_dts_too=False,
                  verbose=False):

    if not hasattr(dts,'__iter__'):
        dts = pd.DatetimeIndex([dts])

    elif not hasattr(dts,'hour'):
        dts = pd.DatetimeIndex(dts)

    if not hasattr(glons,'__iter__'):
        glons = np.array(glons)

    sslat, sslon = map(np.ravel, subsol(dts))
    LTs = ((glon - sslon + 180) % 360 - 180) / 15

    midnightlongitude = (sslon - 180.) % 360
    if return_dts_too:
        return LTs, dts, midnightlongitude, glons
    else:
        return LTs

def cheap_LT_calc(dts,glons,
                  return_dts_too=False,
                  verbose=False):

    # 2021/01/28 NOW we just use Kalle's way of doing things, which uses subsolar longitude
    # if return_dts_too, you get (LTs, dts, midnightlongitude, glons)
    # if not, you get LTs
    print("Using Kalle's subsolar point calculation instead of your junk LT calc!")
    outs = get_localtime(dts,glons,
                         return_dts_too=return_dts_too,
                         verbose=verbose)
    return outs

    # if not hasattr(dts,'__iter__'):
    #     dts = pd.DatetimeIndex([dts])

    # elif not hasattr(dts,'hour'):
    #     dts = pd.DatetimeIndex(dts)

    # if not hasattr(glons,'__iter__'):
    #     glons = np.array(glons)

    # glons = (glons+360) % 360

    # if verbose:
    #     if glons.size ==1:
    #         # relstr = "ahead of"
    #         reltogreenwich = glons/15.

    #         relstr = "ahead of" if (glons <= 180) else "behind"
    #         if reltogreenwich > 12:
    #             reltogreenwich -= 24

    #         print("Longitude {:.2f} is {:.2f} hours {:s} Greenwich!".format(glons,reltogreenwich,relstr))

    # midnightlongitude = -15*(dts.hour.values+dts.minute.values/60+dts.second.values/3600.)
    # midnightlongitude = (midnightlongitude + 360) % 360

    # LTs = (((glons-midnightlongitude) + 360) % 360)/15
    # if return_dts_too:
    #     return LTs, dts, midnightlongitude, glons
    # else:
    #     return LTs

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

    print("Using Kalle's superior subsolar longitude calc")
    sslat, sslon = map(np.ravel, subsol(dts))
    return sslon


    # if not hasattr(dts,'__iter__'):
    #     dts = pd.DatetimeIndex([dts])

    # elif not hasattr(dts,'hour'):
    #     dts = pd.DatetimeIndex(dts)

    # fracHour = dts.hour.values+dts.minute.values/60+dts.second.values/3600.

    # assert not any((fracHour < 0) | (fracHour > 24))

    # fracHour[fracHour > 12] -= 24

    # fracHour *= 15

    # if verbose:
    #     print("Min fracHour: {:.2f}".format(np.min(fracHour)))
    #     print("Max fracHour: {:.2f}".format(np.max(fracHour)))

    # return 180 - fracHour 

