import numpy as np
from hatch_python_utils.earth.geodesy import geo_dist

def crosstrackdistance(lon1,lat1,
                       lon2,lat2,
                       lon3,lat3,
                       rearth=6.370949e3,
                       altitude=0,
                       deg=True,
                       DEBUG=False):
    """
    Distance between the great circle defined by two points P1 and P2
    and a third point P3

    TEST
    ====
    lon1,lat1 = 0,0
    lon2,lat2 = 10,0
    lon3,lat3 = 90,0
    # Distance should be zero
    dist = crosstrackdistance(lon1,lat1,
                       lon2,lat2,
                       lon3,lat3,
                       DEBUG=True)

    REFERENCE
    ==============
    Reproduced from Chris Veness' lovely lon/lat tools via MIT license
    http://www.movable-type.co.uk/scripts/latlong.html
    """

    rearth += altitude

    if deg:
        lat1, lat2, lat3, lon1, lon2, lon3 = map(np.deg2rad, (lat1, lat2, lat3,
                                                              lon1, lon2, lon3))
    
    # Angular distance between point 1 and point 3 
    d13 = geo_dist(lon1,lat1,lon3,lat3,
                   deg=False,
                   rearth=rearth,
                   altitude=0)/rearth

    #initial bearing, 1 to 2 (in DEGREES)
    b12 = bearing(lat1,lat2,lon1,lon2,deg=False)

    #initial bearing, 1 to 3 (in DEGREES)
    b13 = bearing(lat1,lat3,lon1,lon3,deg=False)

    b12,b13 = map(np.deg2rad, (b12,b13))

    dist = np.arcsin(np.sin(d13)*np.sin(b13-b12))*rearth

    if DEBUG:
        print("lon1, lat1   : {:5.2f}°, {:5.2f}°".format(np.rad2deg(lon1),np.rad2deg(lat1)))
        print("lon2, lat2   : {:5.2f}°, {:5.2f}°".format(np.rad2deg(lon2),np.rad2deg(lat2)))
        print("lon3, lat3   : {:5.2f}°, {:5.2f}°".format(np.rad2deg(lon3),np.rad2deg(lat3)))
        print("Sphere radius: {:7.2g}".format(rearth))

        print("Distance between P3 and great circle formed by P1 and P2: ")
        print("{:7.2g}".format(dist))

        print("Angular distance (deg) between P3 and great circle formed by P1 and P2: ")
        print("{:7.2g}".format(np.rad2deg(dist/rearth)))

    return dist


def bearing(lat1, lat2, lon1, lon2, deg=True):
    if deg:
        lat1, lat2, lon1, lon2 = map(np.deg2rad, (lat1, lat2, lon1, lon2))
    x = np.cos(lat2)*np.sin(lon2-lon1)
    y = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    return np.rad2deg(np.arctan2(x, y))


def mlt_bearing(mlats, mlts, shift_to_90deg=False):

    dlat = np.diff(mlats, append=mlats[-1]+(mlats[-1]-mlats[-2]))

    lonunwrap = np.rad2deg(np.unwrap(np.deg2rad(mlts*15.)))
    dlon = np.diff(
        lonunwrap, append=lonunwrap[-1]+(lonunwrap[-1]-lonunwrap[-2]))

    # FEIL
    # np.diff(mlts).max()
    # > 23.999415079752602
    # RIKTIG
    # np.diff(np.rad2deg(lonunwrap)/15.).max()
    # >0.27708638509115247

    lat1 = mlats
    lat2 = mlats+dlat

    lon1 = lonunwrap
    lon2 = lonunwrap+dlon

    if shift_to_90deg:
        magbearing = np.abs(bearing(lat1, lat2, lon1, lon2, deg=True))
        magbearing[magbearing > 90] -= 180
        magbearing = np.abs(magbearing)
        return magbearing
        # mag.loc[(mag['magbearingShift'] > 90),'magbearingshift'] = mag[mag['magbearingShift'] > 90]['magbearingShift']-180
        # mag.loc[mag['magbearingShift'] > 90,'magbearingShift'] -= 180
    else:
        return bearing(lat1, lat2, lon1, lon2, deg=True)
