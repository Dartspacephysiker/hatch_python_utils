# 2019/05/30
from datetime import datetime
import apexpy
import igrf12
from pytt.earth import geodesy
import numpy as np
import pandas as pd


def babyFunc2(a, mlon, datime):
    if np.isnan(mlon):
        mlt = np.nan
    else:
        mlt = a.mlon2mlt(mlon, datime)
    return mlt


def geodetic2apex(*args,
                  apexRefTime=datetime(2012, 1, 1),
                  apexRefHeight_km=110,
                  quiet=False,
                  nancheck=False,
                  returnPandas=True):
    """
    geodetic2apex(gdlat, gdlon, gdalt_km[, times])
    gdlat, gdlon in degrees
    times: datetime list object
    apexRefTime: datetime object
    """

    # assert len(args) >= 3, "geodetic2apex(gdlat, gdlon, gdalt_km[,times])"
    assert len(args) == 4, "geodetic2apex(gdlat, gdlon, gdalt_km,times)"
    gdlat = args[0]
    gdlon = args[1]
    gdalt_km = args[2]

    # print(gdlat)
    # print(gdlon)
    # print(gdalt_km)

    canDoMLT = False
    if len(args) > 3:
        times = args[3]
        canDoMLT = True

    a = apexpy.Apex(apexRefTime, refh=apexRefHeight_km)

    mlat, mlon = a.geo2apex(
        gdlat, gdlon, gdalt_km)
    # This can be replaced with PyAMPS code
    # /SPENCEdata/Research/Satellites/Swarm/pyAMPS/pyamps/mlt_utils.py

    if canDoMLT:
        mlt = np.array([babyFunc2(a, mlonna, datime) for mlonna, datime in zip(
            mlon, times)])

    # quiet = False
    isv = 0
    itype = 1
    if not quiet:
        print("Getting IGRF ...")

    # glat, glon: geographic Latitude, Longitude
    # HAD TO KLUGE gridigrf12 TO GET IT TO WORK: ORIG ONLY USED ONE ALTITUDE

    igrf = igrf12.gridigrf12(times,
                             glat=geodesy.geodetic2geocentriclat(gdlat),
                             glon=gdlon,
                             alt_km=gdalt_km, isv=isv, itype=itype)

    igrfMag = np.sqrt(igrf.east.values**2.+igrf.north.values **
                      2. + igrf.down.values**2.)

    gdlatJ, altJ, XIGRF, ZIGRF = geodesy.geoc2geod(90.-geodesy.geodetic2geocentriclat(gdlat),
                                                   geodesy.geodeticheight2geocentricR(
                                                       gdlat, gdalt_km),
                                                   -igrf.north.values, -igrf.down.values)

    if nancheck:
        nanners = np.isnan(gdlat) | np.isnan(gdlon) | np.isnan(gdalt_km)
        if nanners[nanners].size/nanners.size > 0.0:
            if not quiet:
                print("nannies!"+"{0}".format(nanners[nanners].size))
            gdlat[nanners] = 0
            gdalt_km[nanners] = 0

    if not quiet:
        print("Getting Apex basevectors ...")
    # From Laundal and Richmond (2016):
    # e1 "points eastward along contours of constant λma,"
    # e2 "points equatorward along contours of constant φ ma (magnetic meridians)"
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = a.basevectors_apex(
        gdlat, gdlon, gdalt_km, coords='geo')

    mapratio = 1. / np.linalg.norm(np.cross(d1.T, d2.T), axis=1)

    # del f1, f2, f3, g1, g2, g3, d1, d2, d3

    if nancheck:
        if nanners[nanners].size/nanners.size > 0.0:
            gdlat[nanners] = np.nan
            gdalt_km[nanners] = np.nan

    # dfMInterp.loc[:, 'e10'] = e1[0, ]
    # dfMInterp.loc[:, 'e11'] = e1[1, ]
    # dfMInterp.loc[:, 'e12'] = e1[2, ]
    # dfMInterp.loc[:, 'e20'] = e2[0, ]
    # dfMInterp.loc[:, 'e21'] = e2[1, ]
    # dfMInterp.loc[:, 'e22'] = e2[2, ]

    # output?
    # mlat,mlon,mapratio[,mlt],igrf.east,XIGRF(northGeodetic),-ZIGRF(upGeodetic),igrfMag
    returnList = [mlat, mlon, mapratio]
    rListNames = ['mlat', 'mlon', 'mapratio']

    if canDoMLT:
        returnList.append(mlt)
        rListNames.append('mlt')

    returnList = returnList + [igrf.east.values, XIGRF, -ZIGRF, igrfMag]
    rListNames = rListNames + ['igrfE', 'igrfN', 'igrfU', 'igrfMag']

    if returnPandas:
        df = pd.DataFrame(data=np.vstack(returnList).T, columns=rListNames)
        if canDoMLT:
            df.set_index(times, inplace=True)
        return df

    else:
        return returnList, rListNames
