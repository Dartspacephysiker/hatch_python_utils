# 2019/05/30
from datetime import datetime
import apexpy
import igrf12
from pytt.earth import geodesy
import numpy as np
import pandas as pd
import scipy.io as sio
from hatch_python_utils import arrays as hArr
from hatch_python_utils.pandas_utils import interp_over_nans
from dateutil.relativedelta import relativedelta
from hatch_python_utils.date_time import toYearFraction


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
                  max_N_months_twixt_apexRefTime_and_obs=3,
                  min_time_resolution__sec=1,
                  interpolateArgs={'method': 'time', 'limit': 21},
                  returnPandas=True,
                  return_IGRF=False,
                  return_apex_d_basevecs=False,
                  return_apex_e_basevecs=False,
                  return_apex_f_basevecs=False,
                  return_apex_g_basevecs=False,
                  return_mapratio=False):
    """
    geodetic2apex(gdlat, gdlon, gdalt_km[, times])
    gdlat, gdlon in degrees
    times: datetime list object
    apexRefTime: datetime object
    Returns
    """

    get_apex_basevecs = return_apex_d_basevecs or return_apex_e_basevecs or \
        return_apex_f_basevecs or return_apex_g_basevecs or return_mapratio

    if max_N_months_twixt_apexRefTime_and_obs is None:
        max_N_months_twixt_apexRefTime_and_obs = 0

    # assert len(args) >= 3, "geodetic2apex(gdlat, gdlon, gdalt_km[,times])"
    assert len(args) == 4, "geodetic2apex(gdlat, gdlon, gdalt_km,times)"

    # gdlat = args[0]
    # gdlon = args[1]
    # gdalt_km = args[2]

    # canDoMLT = False
    # if len(args) > 3:
    #     times = args[3]
    #     canDoMLT = True

    canDoMLT = True

    # df = pd.DataFrame({'gdlat':gdlat,'gdlon':gdlon,'gdalt_km':gdalt_km},index=times)

    df = pd.DataFrame(
        {'gdlat': args[0], 'gdlon': args[1], 'gdalt_km': args[2]}, index=args[3])

    # Interp over nans
    checkCols = ['gdlat', 'gdalt_km']
    interp_over_nans(df, checkCols,
                     max_Nsec_twixt_nans=1,
                     max_Nsec_tot=5,
                     interpolateArgs={'method': 'time', 'limit': 21})

    if df.isna().any().any() and not nancheck:
        print("HELP (or just set nancheck=True)!")
        breakpoint()

    ########################################
    # check if we need to downsample/consider subset
    debug = True

    period_df = ((df.iloc[1].name -
                  df.iloc[0].name)/pd.Timedelta('1s'))

    is_subsampled = False
    if period_df < min_time_resolution__sec:

        strider = int(min_time_resolution__sec/period_df)

        if debug:
            print("DEBUG   NEDSAMPLE: Reduserer antall konversjoner med en faktor på {:d}".format(
                strider))

        dfSub = df.iloc[0::strider]
        is_subsampled = True

    else:
        dfSub = df

    # Return these to NaNs further down ...
    if nancheck:

        nanners = dfSub['gdlat'].isna(
        ) | dfSub['gdlon'].isna() | dfSub['gdalt_km'].isna()

        if nanners[nanners].size/nanners.size > 0.0:
            if not quiet:
                print("nannies!"+"{0}".format(nanners[nanners].size))
            dfSub.loc[nanners, 'gdlat'] = 0
            dfSub.loc[nanners, 'gdalt_km'] = 0

    a = apexpy.Apex(apexRefTime, refh=apexRefHeight_km)

    mlat, mlon = a.geo2apex(
        dfSub['gdlat'].values, dfSub['gdlon'].values, dfSub['gdalt_km'].values)
    # This can be replaced with PyAMPS code
    # /SPENCEdata/Research/Satellites/Swarm/pyAMPS/pyamps/mlt_utils.py

    if canDoMLT:

        times = dfSub.index.to_pydatetime()

        if max_N_months_twixt_apexRefTime_and_obs == 0:

            mlt = np.array([babyFunc2(a, mlonna, datime)
                            for mlonna, datime in zip(mlon,
                                                      times)])

        else:

            print("Updating apexRefTime as we go ...")

            mlt = np.zeros(mlat.shape)*np.nan

            wasSorted = hArr.isSorted_alt(times)
            if not wasSorted:
                print("Times not sorted! Sorting ...")
                sortinds = np.argsort(notsorted)
                unsortinds = np.argsort(sortinds)

                times = np.array(times)[sortinds]
                mlat = mlat[sortinds]
                mlon = mlon[sortinds]

            # Now group by max_N_months_twixt_apexRefTime_and_obs
            apexRefTime = times[0]

            relDelta = relativedelta(
                months=max_N_months_twixt_apexRefTime_and_obs)

            maxIterHere = 3000
            nIter = 0
            while apexRefTime < times[-1]:

                # See if we have any here; if not, skip

                ind_timesHere = (times >= apexRefTime) & (
                    times < (apexRefTime+relDelta))

                nIndsHere = np.where(ind_timesHere)[0].size

                if debug:
                    print("DEBUG   {:s} to {:s} : Got {:d} inds for MLT conversion".format(
                        apexRefTime.strftime("%Y%m%d"),
                        (apexRefTime+relDelta).strftime("%Y%m%d"),
                        nIndsHere))

                if nIndsHere == 0:
                    # Increment apexRefTime by relDelta
                    apexRefTime += relDelta
                    continue

                a.set_epoch(toYearFraction(apexRefTime))

                maxNIndsSamtidig = 500000
                if nIndsHere > maxNIndsSamtidig:
                    print("Break it up...")

                    indbatchCounter = 0
                    nIndsConverted = 0
                    while nIndsConverted < nIndsHere:

                        startIndInd = nIndsConverted
                        stopIndInd = np.min(
                            [startIndInd+maxNIndsSamtidig, nIndsHere])
                        nToConvert = stopIndInd-startIndInd

                        tmpUseInds = np.where(ind_timesHere)[
                            0][startIndInd:stopIndInd]
                        mlt[tmpUseInds] = np.array([babyFunc2(a, mlonna, datime)
                                                    for mlonna, datime in zip(mlon[tmpUseInds],
                                                                              times[tmpUseInds])])

                        nIndsConverted += nToConvert

                else:
                    mlt[ind_timesHere] = np.array([babyFunc2(a, mlonna, datime)
                                                   for mlonna, datime in zip(mlon[ind_timesHere],
                                                                             times[ind_timesHere])])

                # Increment apexRefTime by relDelta
                apexRefTime += relDelta

                nIter += 1
                if nIter >= maxIterHere:
                    print("Too many iterations! Breaking ...")
                    break

            if not wasSorted:
                print("Unsorting things again, you filthy animal")
                times = list(times[unsortinds])
                mlat = mlat[unsortinds]
                mlon = mlon[unsortinds]
                mlt = mlt[unsortinds]

    # returnList = [mlat, mlon, mapratio]
    # rListNames = ['mlat', 'mlon', 'mapratio']

    returnList = [mlat, mlon]
    rListNames = ['mlat', 'mlon']

    if canDoMLT:
        returnList.append(mlt)
        rListNames.append('mlt')

    if return_IGRF:
        # quiet = False
        isv = 0
        itype = 1
        if not quiet:
            print("Getting IGRF ...")

        # glat, glon: geographic Latitude, Longitude
        # HAD TO KLUGE gridigrf12 TO GET IT TO WORK: ORIG ONLY USED ONE ALTITUDE

        igrf = igrf12.gridigrf12(times,
                                 glat=geodesy.geodetic2geocentriclat(
                                     dfSub['gdlat'].values),
                                 glon=dfSub['gdlon'].values,
                                 alt_km=dfSub['gdalt_km'].values, isv=isv, itype=itype)

        igrfMag = np.sqrt(igrf.east.values**2.+igrf.north.values **
                          2. + igrf.down.values**2.)

        gdlatJ, altJ, XIGRF, ZIGRF = geodesy.geoc2geod(90.-geodesy.geodetic2geocentriclat(dfSub['gdlat'].values),
                                                       geodesy.geodeticheight2geocentricR(
                                                           dfSub['gdlat'].values, dfSub['gdalt_km'].values),
                                                       -igrf.north.values, -igrf.down.values)

        # output?
        # mlat,mlon,mapratio[,mlt],igrf.east,XIGRF(northGeodetic),-ZIGRF(upGeodetic),igrfMag

        returnList = returnList + [igrf.east.values, XIGRF, -ZIGRF, igrfMag]
        rListNames = rListNames + ['igrfE', 'igrfN', 'igrfU', 'igrfMag']

    if get_apex_basevecs:

        if not quiet:
            print("Getting Apex basevectors ...")
        # From Laundal and Richmond (2016):
        # e1 "points eastward along contours of constant λma,"
        # e2 "points equatorward along contours of constant φ ma (magnetic meridians)"
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = a.basevectors_apex(
            dfSub['gdlat'].values, dfSub['gdlon'].values, dfSub['gdalt_km'].values, coords='geo')

        if return_mapratio:
            mapratio = 1. / np.linalg.norm(np.cross(d1.T, d2.T), axis=1)
            returnList = returnList + [mapratio]
            rListNames = rListNames + ['mapratio']

        # del f1, f2, f3, g1, g2, g3, d1, d2, d3

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

    if nancheck:
        if nanners[nanners].size/nanners.size > 0.0:
            dfSub.loc[nanners, 'gdlat'] = np.nan
            dfSub.loc[nanners, 'gdalt_km'] = np.nan

    ########################################
    # Make final outputdict
    returnDict = {key: val for key, val in zip(rListNames, returnList)}

    if returnPandas:

        if is_subsampled:
            intoApexIndex = df.index.intersection(dfSub.index)

            dfOut = pd.DataFrame(columns=rListNames, dtype=np.float64,
                                 index=df.index)

            # dfOut = pd.DataFrame(returnDict,
            #                      index=df.index)

            # dfOut['utc'] = dfOut.index
            # dfOut.index = (dfOut.index-dfOut.index[0])/pd.Timedelta('1s')
            # interpolateArgs['method'] = 'index'
            # if 'limit' in interpolateArgs.keys():
            #     _ = interpolateArgs.pop('limit')

            shouldBeUnwrapped = ['mlon', 'mlt']
            for col in rListNames:

                if col in shouldBeUnwrapped:
                    if col == 'mlon':
                        dfOut.loc[intoApexIndex, col] = np.unwrap(
                            np.deg2rad(returnDict['mlon']))
                    elif col == 'mlt':
                        dfOut.loc[intoApexIndex, 'mlt'] = np.unwrap(
                            np.deg2rad(returnDict['mlt']*15.))
                    else:
                        print("What is this?")
                        breakpoint()

                else:
                    dfOut.loc[intoApexIndex, col] = returnDict[col]

                dfOut.loc[:, col] = dfOut[col].interpolate(**interpolateArgs)

                if col in shouldBeUnwrapped:
                    if col == 'mlon':
                        dfOut.loc[:, 'mlon'] = np.rad2deg(
                            (dfOut['mlon'].values + np.pi) % (2 * np.pi) - np.pi)
                    elif col == 'mlt':
                        dfOut.loc[:, 'mlt'] = np.rad2deg(
                            (dfOut['mlt'].values) % (2 * np.pi))/15.

        else:
            dfOut = pd.DataFrame(data=np.vstack(returnList).T,
                                 columns=rListNames, index=df.index)

        return dfOut

    else:
        return returnDict


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

        tmpea = sio.readsav(self._indir+self._infile, python_dict=True)['ea']

        newea = np.zeros(tmpea[0][0].shape, np.dtype(dict(names=['mini', 'maxi', 'centeri', 'minm', 'maxm', 'centerm'],
                                                          formats=[np.float32, np.float32, np.float32,
                                                                   np.float32, np.float32, np.float32])))
        newea['mini'] = tmpea['mini'].item()
        newea['maxi'] = tmpea['maxi'].item()
        newea['centeri'] = (tmpea.maxi.item()+tmpea.mini.item())/2.
        newea['minm'] = tmpea['minm'].item()
        newea['maxm'] = tmpea['maxm'].item()
        newea['centerm'] = (tmpea.maxm.item()+tmpea.minm.item())/2.

        self.ea = newea.view(np.recarray)

        self._hemi == hemi

        if hemi.lower() == 'south':

            print("Using south stuff")
            # tmpMinI = (-1.)*np.flip(self.ea.maxi)
            # tmpMaxI = (-1.)*np.flip(self.ea.mini)

            # tmpMinM = np.flip(self.ea.minm)
            # tmpMaxM = np.flip(self.ea.maxm)

            tmpMinI = (-1.)*self.ea.maxi
            tmpMaxI = (-1.)*self.ea.mini

            tmpMinM = self.ea.minm
            tmpMaxM = self.ea.maxm

            self.ea.mini = tmpMinI
            self.ea.maxi = tmpMaxI
            self.ea.centeri = (tmpMinI+tmpMaxI)/2.

            self.ea.minm = tmpMinM
            self.ea.maxm = tmpMaxM
            self.ea.centerm = (tmpMinM+tmpMaxM)/2.


def geoclatR2geodlatheight(glat, r_km):
    """
    glat : geocentric latitude (degrees)
    r_km : radius (km)
    Returns gdlat, gdalt_km
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
    gdalt_km = r_km * np.cos(DLTCL) - a * np.sqrt(1 - E2 * np.sin(gdlat) ** 2)

    gdlat = gdlat / d2r

    return gdlat, gdalt_km


def bin_into_equal_area(lats, mlts, data,
                        statfunc=np.median,
                        ea=None,
                        verbose=False,
                        add_count=False,
                        return_indices=False):
    """
    'data' is assumed to have rows of observations, and columns of different observation types
    """

    isPandas = isinstance(data, pd.core.frame.DataFrame) or \
        isinstance(data, pd.core.series.Series)

    if isPandas:
        def intrastatfunc(data, inds):
            return statfunc(data[inds], axis=0)
    else:
        def intrastatfunc(data, inds):
            return statfunc(data[inds, :], axis=0)

    if ea is None:
        print("Loading in equal-area thing sjølv")
        ea = EqualAreaBins()

    latII = pd.IntervalIndex.from_arrays(ea.ea.mini, ea.ea.maxi, closed='left')
    mltII = pd.IntervalIndex.from_arrays(ea.ea.minm, ea.ea.maxm, closed='left')

    # indlist = []
    statlist = []

    if return_indices:
        returners = []

    assert len(data.shape) < 3, "Can't handle three-dim data!"

    if len(data.shape) == 2:
        nObs = data.shape[0]
        nCols = data.shape[1]

        # if add_count:
        #     nCols += 1

        if (nCols > 20) or ((nCols/nObs) > 10):
            print("Too funky, not sure what to do with your data")
            breakpoint()

    else:
        nObs = data.shape[0]
        nCols = 1

    # if add_count:
    #     nCols += 1

    if add_count:
        nCols += 1
        addcountcol = nCols-1
        lastdataind = nCols-1
    else:
        lastdataind = nCols

    if nCols > 1:
        blankrow = np.array([np.nan]*nCols)
    else:
        blankrow = np.nan

    statlist = np.zeros((ea.ea.maxi.size, nCols), dtype=np.float64)*np.nan

    for i, (latbin, mltbin) in enumerate(zip(latII, mltII)):
        # tmpstat = blankrow.copy()
        if verbose:
            print(i, latbin, mltbin)

        inds = (lats >= latbin.left) \
            & (lats < latbin.right) \
            & (mlts >= mltbin.left) \
            & (mlts < mltbin.right)
        # indlist.append(np.where(inds)[0])
        nHere = np.where(inds)[0].size
        if nHere > 0:
            # breakpoint()
            # tmpstat = statfunc(data[inds,:],axis=0)
            # tmpstat[:nCols-add_count] = intrastatfunc(data, inds)
            statlist[i, :lastdataind] = intrastatfunc(data, inds)

        if add_count:
            # tmpstat[-1] = nHere
            # print(i, latbin, mltbin, nHere)
            statlist[i, addcountcol] = nHere

        if return_indices:
            returners.append(np.where(inds)[0])

        # statlist.append(tmpstat)

    # statlist = np.array(statlist)

    if return_indices:
        return statlist, return_indices
    else:
        return statlist


def get_magnetic_polcap_equalarea_bin_weights(ea, apexObj, polcaplowlat=70.,
                                              ea_altitude=0,
                                              mirror_SH=False,
                                              verbose=False):

    isSH = ea.hemi.lower() == 'south'
    if isSH and not mirror_SH:
        if polcaplowlat > 0:
            polcaplowlat *= -1.
        # print("Is southern, wif polcaplowlat {:.2f}".format(polcaplowlat))

        def tmpcompfunc(ea_mlat, polcaplowlat):
            return ea_mlat <= polcaplowlat
    else:
        def tmpcompfunc(ea_mlat, polcaplowlat):
            return ea_mlat >= polcaplowlat

    nboxes = 9
    ea_alts = np.array([ea_altitude]*len(ea.ea.maxi))

    # weights = np.zeros((ea.ea.maxi.size,nboxes))
    weights = np.zeros(ea.ea.maxi.size, dtype=np.float64)

    # Divide each into 9 boxes, count up how many
    # LON positions: LEFT, LEFT-CENTER, CENTER, CENTER-RIGHT, RIGHT [L, LC, C, CR, R]
    # LAT positions: BOTTOM, BOTTOM-MID, MID, MID-TOP, TOP          [B, BM, M, MT, T]

    # 'ea_9sq' = "Equal-area divided into nine squares"
    ea_9sq = pd.DataFrame(dict(latbm=(ea.ea.mini+ea.ea.centeri)/2,
                               latm=ea.ea.centeri,
                               latmt=(ea.ea.maxi+ea.ea.centeri)/2,
                               mltlc=(ea.ea.minm+ea.ea.centerm)/2,
                               mltc=ea.ea.centerm,
                               mltcr=(ea.ea.maxm+ea.ea.centerm)/2))

    # Top row
    box0 = [ea_9sq['mltlc'].values*15., ea_9sq['latmt'].values]
    box1 = [ea_9sq['mltc'].values*15., ea_9sq['latmt'].values]
    box2 = [ea_9sq['mltcr'].values*15., ea_9sq['latmt'].values]

    # Mid row
    box3 = [ea_9sq['mltlc'].values*15., ea_9sq['latm'].values]
    box4 = [ea_9sq['mltc'].values*15., ea_9sq['latm'].values]
    box5 = [ea_9sq['mltcr'].values*15., ea_9sq['latm'].values]

    # Bottom row
    box6 = [ea_9sq['mltlc'].values*15., ea_9sq['latbm'].values]
    box7 = [ea_9sq['mltc'].values*15., ea_9sq['latbm'].values]
    box8 = [ea_9sq['mltcr'].values*15., ea_9sq['latbm'].values]

    # breakpoint()

    for ibox, box in enumerate([box0, box1, box2,
                                box3, box4, box5,
                                box6, box7, box8]):
        if verbose:
            print("Box #{:d}".format(ibox))

        ea_mlat, ea_mlon = apexObj.geo2apex(box[1], box[0], ea_alts)
        # incap_inds = ea_mlat >= polcaplowlat
        incap_inds = tmpcompfunc(ea_mlat, polcaplowlat)
        weights[incap_inds] += 1.

    weights = weights/nboxes

    return weights
