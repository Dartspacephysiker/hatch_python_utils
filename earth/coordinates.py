# 2019/05/30
from datetime import datetime
import apexpy
import igrf12
from pytt.earth import geodesy
import numpy as np
import pandas as pd
import scipy.io as sio


def babyFunc2(a, mlon, datime):
    if np.isnan(mlon):
        mlt = np.nan
    else:
        mlt = a.mlon2mlt(mlon, datime)
    return mlt


def ECEFtoENUMatrix(lon, gdlat):

    lamb = np.deg2rad(lon)
    phi = np.deg2rad(gdlat)

    return np.stack([np.vstack([-np.sin(lamb), np.cos(lamb), np.zeros(lamb.size)]),
                     np.vstack([-np.sin(phi)*np.cos(lamb), -
                                np.sin(phi)*np.sin(lamb), np.cos(phi)]),
                     np.vstack([np.cos(phi)*np.cos(lamb), np.cos(phi)*np.sin(lamb), np.sin(phi)])])

    # return np.vstack([[-np.sin(lamb)            , np.cos(lamb)            ,0          ],
    #                   [-np.sin(phi)*np.cos(lamb),-np.sin(phi)*np.sin(lamb),np.cos(phi)],
    #                   [ np.cos(phi)*np.cos(lamb), np.cos(phi)*np.sin(lamb),np.sin(phi)]])


def geodetic2apex(*args,
                  apexRefTime=datetime(2012, 1, 1),
                  apexRefHeight_km=110,
                  quiet=False,
                  nancheck=False,
                  returnPandas=True,
                  return_apex_d_basevecs=False,
                  return_apex_e_basevecs=False,
                  return_apex_f_basevecs=False,
                  return_apex_g_basevecs=False):
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

    returnList = [mlat, mlon, mapratio]
    rListNames = ['mlat', 'mlon', 'mapratio']

    # if return_apex_basevecs:

    if return_apex_d_basevecs:
        returnList = returnList + [d1[0, ], d1[1, ], d1[2, ],
                                   d2[0, ], d2[1, ], d2[2, ],
                                   d3[0, ], d3[1, ], d3[2, ]]
        rListNames = rListNames + ['d10', 'd11', 'd12',
                                   'd20', 'd21', 'd22',
                                   'd30', 'd31', 'd32']

    if return_apex_e_basevecs:
        returnList = returnList + [e1[0, ], e1[1, ], e1[2, ],
                                   e2[0, ], e2[1, ], e2[2, ],
                                   e3[0, ], e3[1, ], e3[2, ]]
        rListNames = rListNames + ['e10', 'e11', 'e12',
                                   'e20', 'e21', 'e22',
                                   'e30', 'e31', 'e32']
    if return_apex_f_basevecs:
        returnList = returnList + [f1[0, ], f1[1, ], f1[2, ],
                                   f2[0, ], f2[1, ], f2[2, ],
                                   f3[0, ], f3[1, ], f3[2, ]]
        rListNames = rListNames + ['f10', 'f11', 'f12',
                                   'f20', 'f21', 'f22',
                                   'f30', 'f31', 'f32']
    if return_apex_g_basevecs:
        returnList = returnList + [g1[0, ], g1[1, ], g1[2, ],
                                   g2[0, ], g2[1, ], g2[2, ],
                                   g3[0, ], g3[1, ], g3[2, ]]
        rListNames = rListNames + ['g10', 'g11', 'g12',
                                   'g20', 'g21', 'g22',
                                   'g30', 'g31', 'g32']

    # dfMInterp.loc[:, 'e10'] = e1[0, ]
    # dfMInterp.loc[:, 'e11'] = e1[1, ]
    # dfMInterp.loc[:, 'e12'] = e1[2, ]
    # dfMInterp.loc[:, 'e20'] = e2[0, ]
    # dfMInterp.loc[:, 'e21'] = e2[1, ]
    # dfMInterp.loc[:, 'e22'] = e2[2, ]

    # output?
    # mlat,mlon,mapratio[,mlt],igrf.east,XIGRF(northGeodetic),-ZIGRF(upGeodetic),igrfMag

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


class EqualAreaBins(object):
    """
    See example usage in journal__20190815__lek_med_Swarm_crosstrack.ipynb
    """

    def __init__(self, hemi='north'):

        self._indir = '/SPENCEdata/Research/database/equal-area_binning/'
        self._infile = 'equalArea--20161014--struct_and_ASCII_tmplt.idl'

        self.ea = None
        self._hemi = hemi

        self.set_hemi(hemi)

        # self.ea = sio.readsav(indir+infile, python_dict=True)['ea']

        # self._hemi = hemi

    @property
    def hemi(self):
        """Get the current hemi."""
        return self._hemi

    def set_hemi(self, hemi):

        allowedHemis = ['north', 'south']

        if hemi not in allowedHemis:
            print("Must select either 'north' or 'south'!")
            return

        self.ea = sio.readsav(self._indir+self._infile, python_dict=True)['ea']

        self._hemi == hemi

        if hemi.lower() == 'south':

            print("Using south stuff (never tested!)")
            tmpMinI = (-1.)*np.flip(self.ea.maxi)
            tmpMaxI = (-1.)*np.flip(self.ea.mini)

            tmpMinM = np.flip(self.ea.minm)
            tmpMaxM = np.flip(self.ea.maxm)

            self.ea.mini = tmpMinI
            self.ea.maxi = tmpMaxI

            self.ea.minm = tmpMinM
            self.ea.maxm = tmpMaxM


def geoclatR2geodlatheight(glat, r_km):
    """
    glat : geocentric latitude (degrees)
    r_km : radius (km)
    Ripped off Kalle's pytt.geodesy.geoc2geod ...
    """
    d2r = np.pi/180
    r2d = 180 / np.pi
    WGS84_e2 = 0.00669437999014
    WGS84_a = 6378.137

    a = WGS84_a
    b = a*np.sqrt(1 - WGS84_e2)

    # breakpoint()

    E2 = 1.-(b/a)**2
    E4 = E2*E2
    E6 = E4*E2
    E8 = E4*E4
    A21 = (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024.
    A22 = (E6 + E8) / 32.
    A23 = -3.*(4.*E6 + 3.*E8) / 256.
    A41 = -(64.*E4 + 48.*E6 + 35.*E8)/1024.
    A42 = (4.*E4 + 2.*E6 + E8) / 16.
    A43 = 15.*E8 / 256.
    A44 = -E8 / 16.
    A61 = 3.*(4.*E6 + 5.*E8)/1024.
    A62 = -3.*(E6 + E8) / 32.
    A63 = 35.*(4.*E6 + 3.*E8) / 768.
    A81 = -5.*E8 / 2048.
    A82 = 64.*E8 / 2048.
    A83 = -252.*E8 / 2048.
    A84 = 320.*E8 / 2048.

    SCL = np.sin(glat * d2r)

    RI = a/r_km
    A2 = RI*(A21 + RI * (A22 + RI * A23))
    A4 = RI*(A41 + RI * (A42 + RI*(A43+RI*A44)))
    A6 = RI*(A61 + RI * (A62 + RI * A63))
    A8 = RI*(A81 + RI * (A82 + RI*(A83+RI*A84)))

    CCL = np.sqrt(1-SCL**2)
    S2CL = 2.*SCL * CCL
    C2CL = 2.*CCL * CCL-1.
    S4CL = 2.*S2CL * C2CL
    C4CL = 2.*C2CL * C2CL-1.
    S8CL = 2.*S4CL * C4CL
    S6CL = S2CL * C4CL + C2CL * S4CL

    DLTCL = S2CL * A2 + S4CL * A4 + S6CL * A6 + S8CL * A8
    gdlat = DLTCL + glat * d2r
    height = r_km * np.cos(DLTCL) - a * np.sqrt(1 - E2 * np.sin(gdlat) ** 2)

    gdlat = gdlat / d2r

    return gdlat, height
