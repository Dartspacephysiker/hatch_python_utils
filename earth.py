import numpy as np
import igrf12


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


def sphDist(lat1, mlt1, lat2, mlt2):
    """
    Great-circle distance
    lat1: Latitude of 1st point in degrees
    lat2: Latitude of 2nd point in degrees
    mlt1: Magnetic local time of 1st point in hours (so between 0 and 24)
    mlt2: Magnetic local time of 2nd point in hours
    """
    lat1R = np.deg2rad(lat1)
    mlt1R = np.deg2rad(mlt1*15.)
    lat2R = np.deg2rad(lat2)
    mlt2R = np.deg2rad(mlt2*15.)

    return np.rad2deg(np.arccos(np.sin(lat1R)*np.sin(lat2R)+np.cos(lat1R)*np.cos(lat2R)*np.cos(np.abs(mlt1R-mlt2R))))

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
