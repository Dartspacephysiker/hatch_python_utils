import numpy as np
import igrf12
from datetime import datetime


def earthSunDist(doy):
    """
    https://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year

    This is an approximate expression. Term by term,

    1       : The mean distance between the Earth and the Sun is about one
    astronomical unit.

    0.01672 : This is the eccentricity of the Earth's about about the Sun.

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

# def get_igrf_field(datetimes,geolat,geolon,gdalt_km):

#     quiet = False
#     isv = 0
#     itype = 1
#     if not quiet:
#         print("Getting IGRF ...")

#     # glat, glon: geographic Latitude, Longitude
#     # HAD TO KLUGE gridigrf12 TO GET IT TO WORK: ORIG ONLY USED ONE ALTITUDE
#     igrf = igrf12.gridigrf12(datetimes,
#                              # glat=dfMInterp.gdlat,
#                              glat=geolat,
#                              glon=geolon,
#                              alt_km=gdalt_km, isv=isv, itype=itype)

# # geodeticheight2geocentricR(lat, height)


#     gdlatJ, altJ, XIGRF, ZIGRF = geodesy.geoc2geod(90.-geolat,
#                                                    dfMInterp["Radius"].values/1000.,
#                                                    -igrf.north.values, -igrf.down.values)
#     igrfMag = np.sqrt(igrf.east.values**2.+igrf.north.values **
#                       2. + igrf.down.values**2.)

#     return

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
