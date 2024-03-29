# 2019/05/30
from datetime import datetime
import apexpy
# import igrf12 as igrf
from pysymmetry import geodesy
import numpy as np
import pandas as pd
import scipy.io as sio
from hatch_python_utils import arrays as hArr
from hatch_python_utils.pandas_utils import interp_over_nans
from dateutil.relativedelta import relativedelta
from hatch_python_utils.date_time import toYearFraction
from pyamps.mlt_utils import mlon_to_mlt


# def babyFunc2(a, mlon, datime):
#     if np.isnan(mlon):
#         mlt = np.nan
#     else:
#         mlt = a.mlon2mlt(mlon, datime)
#     return mlt


def convapexvectoENU(vtheta,vphi,vpar,mlt,mlat,dt,refh=130):
    """
    vtheta: Vector component in magnetic north/south direction, with equatorward positive (at least in NH)
    vphi  : Vector component in magnetic east/west direction, with eastward positive
    vpar  : Vector component along magnetic field
    dt    : Reference datetime, used for creating Apex object
    """

    a = apexpy.Apex(dt,refh=refh)
    mlon = a.mlt2mlon(mlt,dt)
    glat, glon, error = a.apex2geo(mlat,mlon,refh)
    

    # Apex basevectors - components are east, north, and up
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = a.basevectors_apex(mlat.flatten(),
                                                                        mlon.flatten(),
                                                                        refh,
                                                                        coords='apex')

    vec = np.vstack([vphi.flatten(),vtheta.flatten(),vpar.flatten()])

    # first = np.array([np.sum(np.array([d1[0,0],d2[0,0],d3[0,0]])*vec[:,0]),
    #                   np.sum(np.array([d1[1,0],d2[1,0],d3[1,0]])*vec[:,0]),
    #                   np.sum(np.array([d1[2,0],d2[2,0],d3[2,0]])*vec[:,0])])
    # array([-0.03262529,  0.02291903,  0.00150332])

    #One matrix
    # np.vstack([d1[:,0],d2[:,0],d3[:,0]])

    APEX2GEO = np.array([d1,d2,d3])
    
    # NOT SURE IF THIS IS NECESSARY
    # APEX2GEO = np.transpose(APEX2GEO,(1,0,2))  # so that APEX2GEO[:,0,0] selects

    vecENU = np.einsum('ij...,i...',APEX2GEO,vec).T

    return vecENU,glat,glon


def ENUtoECEFMatrix(lon, gclat):
    """
    The transpose of what you get from ECEFtoENUMatrix
    """

    return np.stack([np.vstack([             -np.sin(phi), -np.sin(lamb)*np.cos(phi), np.cos(lamb)*np.cos(phi)]),
                     np.vstack([              np.cos(phi), -np.sin(lamb)*np.sin(phi), np.cos(lamb)*np.sin(phi)]),
                     np.vstack([       np.zeros(phi.size),              np.cos(lamb),             np.sin(lamb)])])



def ECEFtoENUMatrix(lon, gclat):
    """

    First  row is "East"  (phi)      vector (   phihat = -sinp xhat + cosp yhat)
    Second row is "North" (90-theta) vector (lambdahat = -sinl cosp xhat - cosl sinp yhat + cosl zhat)
    Third  row is "Up"    (radial)   vector (     rhat =  cosl cosp xhat + cosl sinp yhat + sinl zhat)

    Here phi    ("p") is the azimuthal coordinate in spherical coordinates
         lambda ("l") is the 90°-theta latitude angle


    Naturally, the columns of this matrix give xhat, yhat, and zhat (in that order) in terms of 
    ehat, nhat, and uhat

    So! Multiply a vector with ENU components by this matrix to get back a vector in ECEF coordinates
    """

    phi = np.deg2rad(lon)
    lamb = np.deg2rad(gclat)

    return np.stack([np.vstack([             -np.sin(phi),               np.cos(phi),  np.zeros(phi.size)]),
                     np.vstack([-np.sin(lamb)*np.cos(phi), -np.sin(lamb)*np.sin(phi),        np.cos(lamb)]),
                     np.vstack([ np.cos(lamb)*np.cos(phi),  np.cos(lamb)*np.sin(phi),        np.sin(lamb)])])

    # return np.vstack([[-np.sin(phi)            , np.cos(phi)            ,0          ],
    #                   [-np.sin(lamb)*np.cos(phi),-np.sin(lamb)*np.sin(phi),np.cos(lamb)],
    #                   [ np.cos(lamb)*np.cos(phi), np.cos(lamb)*np.sin(phi),np.sin(lamb)]])


def geodetic2apex(*args,
                  apexRefTime=datetime(2012, 1, 1),
                  apexRefHeight_km=110,
                  quiet=False,
                  nancheck=False,
                  max_N_months_twixt_apexRefTime_and_obs=3,
                  min_time_resolution__sec=1,
                  interpolateArgs={'method': 'time', 'limit': 21},
                  do_qdcoords=False,
                  returnPandas=True,
                  return_IGRF=False,
                  return_apex_d_basevecs=False,
                  return_apex_e_basevecs=False,
                  return_apex_f_basevecs=False,
                  return_apex_g_basevecs=False,
                  return_dipoletilt=False,
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

    # if get_apex_basevecs or return_mapratio or return_IGRF:
    if return_IGRF:
        # assert 2<0,"WARNING! geodetic2apex takes care of changes in time with mlat, mlon, mlt, and basevectors, but not with mapratio or IGRF!"
        # assert 2<0,"WARNING! geodetic2apex takes care of main field temporal evolution for mlat, mlon, mlt, mapratio, and basevectors, but not IGRF!"
        import warnings
        warnings.warn("New version (2021/07/08) of this function, where return_IGRF involves looping over time, is untested!", UserWarning)


    if max_N_months_twixt_apexRefTime_and_obs is None:
        max_N_months_twixt_apexRefTime_and_obs = 0

    # assert len(args) >= 3, "geodetic2apex(gdlat, gdlon, gdalt_km[,times])"
    assert len(args) == 4, "geodetic2apex(gdlat, gdlon, gdalt_km,times)"

    # Find out --- should we even worry about time here?
    haveDTIndex = isinstance(args[3],pd.DatetimeIndex)
    multitime = False
    if haveDTIndex:
        multitime = True
    else:
        if isinstance(args[3],list) or isinstance(args[3],np.ndarray):
            if len(args[3]) == len(args[2]):
                multitime = True
            # else:
            #     multitime = False

    # gdlat = args[0]
    # gdlon = args[1]
    # gdalt_km = args[2]

    # canDoMLT = False
    # if len(args) > 3:
    #     times = args[3]
    #     canDoMLT = True

    canDoMLT = True

    # Set up apex object
    a = apexpy.Apex(apexRefTime, refh=apexRefHeight_km)

    # If not doing multiple times, just do the conversions and get out
    if not multitime:
        if do_qdcoords:
            mlat, mlon = a.geo2qd(
                args[0], args[1], args[2])
        else:
            mlat, mlon = a.geo2apex(
                args[0], args[1], args[2])
    
        returnList = [mlat, mlon]
        rListNames = ['mlat', 'mlon']
    
        if canDoMLT:
            mlt = mlon_to_mlt(mlon, [args[3]]*len(mlon), args[3].year)
    
            returnList.append(mlt)
            rListNames.append('mlt')

        returnDict = {key: val for key, val in zip(rListNames, returnList)}

        if returnPandas:
            dfOut = pd.DataFrame(data=np.vstack(returnList).T,columns=rListNames)

            return dfOut

        else:
            return returnDict

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


    # Now group by max_N_months_twixt_apexRefTime_and_obs
    times = dfSub.index.to_pydatetime()
    apexRefTime = times[0]

    relDelta = relativedelta(
        months=max_N_months_twixt_apexRefTime_and_obs)

    maxIterHere = 3000
    nIter = 0

    NOBS = dfSub['gdlat'].values.shape[0]
    mlat = np.zeros(NOBS)*np.nan
    mlon = np.zeros(NOBS)*np.nan

    if get_apex_basevecs:
        
        if return_apex_d_basevecs:
            # d10,d11,d12 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # d20,d21,d22 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # d30,d31,d32 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            d1,d2,d3 = np.zeros((3,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan
        if return_apex_e_basevecs:
            # e10,e11,e12 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # e20,e21,e22 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # e30,e31,e32 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            e1,e2,e3 = np.zeros((3,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan
        if return_apex_f_basevecs:
            # f10,f11,f12 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # f20,f21,f22 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # f30,f31,f32 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            f1,f2,f3 = np.zeros((2,NOBS))*np.nan,np.zeros((2,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan
        if return_apex_g_basevecs:
            # g10,g11,g12 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # g20,g21,g22 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            # g30,g31,g32 = np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan,np.zeros(NOBS)*np.nan
            g1,g2,g3 = np.zeros((3,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan,np.zeros((3,NOBS))*np.nan

    if return_IGRF:
        
        if not quiet:
            print("Getting IGRF ...")

        #igrf.east.values, XIGRF, -ZIGRF
        EIGRF,NIGRF,UIGRF,igrfMag = np.zeros(NOBS),np.zeros(NOBS),np.zeros(NOBS),np.zeros(NOBS)

    if return_dipoletilt:
        from dipole import dipole_tilt
        from hatch_python_utils.time_tools import datetime_to_yearfrac
        dptilt = np.zeros(NOBS)*np.nan

    while apexRefTime < times[-1]:

        # See if we have any here; if not, skip

        ind_timesHere = (times >= apexRefTime) & (
            times < (apexRefTime+relDelta))

        nIndsHere = np.where(ind_timesHere)[0].size

        if debug:
            print("DEBUG   {:s} to {:s} : Got {:d} inds for mlat/mlon conversion".format(
                apexRefTime.strftime("%Y%m%d"),
                (apexRefTime+relDelta).strftime("%Y%m%d"),
                nIndsHere))

        if nIndsHere == 0:
            # Increment apexRefTime by relDelta
            apexRefTime += relDelta
            nIter += 1
            continue

        a.set_epoch(toYearFraction(apexRefTime))

        if do_qdcoords:
            mlattmp, mlontmp = a.geo2qd(
                dfSub['gdlat'].values[ind_timesHere],
                dfSub['gdlon'].values[ind_timesHere],
                dfSub['gdalt_km'].values[ind_timesHere])
        else:
            mlattmp, mlontmp = a.geo2apex(
                dfSub['gdlat'].values[ind_timesHere],
                dfSub['gdlon'].values[ind_timesHere],
                dfSub['gdalt_km'].values[ind_timesHere])

        mlat[ind_timesHere] = mlattmp
        mlon[ind_timesHere] = mlontmp

        # Basevectors
        if get_apex_basevecs:

            # if not quiet:
            #     print("Getting Apex basevectors ...")
            # From Laundal and Richmond (2016):
            # e1 "points eastward along contours of constant λma,"
            # e2 "points equatorward along contours of constant φ ma (magnetic meridians)"
            # "t" stands for "temporary"
            #
            # From apexpy.apex.basevectors_apex:
            # "vector components are geodetic east, north, and up (only east and north for `f1` and `f2`)"
            f1t, f2t, f3t, g1t, g2t, g3t, d1t, d2t, d3t, e1t, e2t, e3t = a.basevectors_apex(
                dfSub['gdlat'].values[ind_timesHere],
                dfSub['gdlon'].values[ind_timesHere],
                dfSub['gdalt_km'].values[ind_timesHere], coords='geo')

            if return_apex_d_basevecs or return_mapratio:
                d1[:,ind_timesHere] = d1t
                d2[:,ind_timesHere] = d2t
                d3[:,ind_timesHere] = d3t
            if return_apex_e_basevecs:
                e1[:,ind_timesHere] = e1t
                e2[:,ind_timesHere] = e2t
                e3[:,ind_timesHere] = e3t
            if return_apex_f_basevecs:
                f1[:,ind_timesHere] = f1t
                f2[:,ind_timesHere] = f2t
                f3[:,ind_timesHere] = f3t
            if return_apex_g_basevecs:
                g1[:,ind_timesHere] = g1t
                g2[:,ind_timesHere] = g2t
                g3[:,ind_timesHere] = g3t


        if return_IGRF:

            isv = 0
            itype = 1
            
            # glat, glon: geographic Latitude, Longitude
            # HAD TO KLUGE gridigrf12 TO GET IT TO WORK: ORIG ONLY USED ONE ALTITUDE
            
            # OLD CALL TO igrf12 MODULE
            # igrf = igrf12.gridigrf12(times[ind_timesHere],
            #                          glat=geodesy.geodetic2geocentriclat(
            #                              dfSub['gdlat'].values[ind_timesHere]),
            #                          glon=dfSub['gdlon'].values[ind_timesHere],
            #                          alt_km=dfSub['gdalt_km'].values[ind_timesHere], isv=isv, itype=itype)
            
            # NEW CALL??
            # breakpoint()
            # igrfdf = igrf12.gridigrf12(times[ind_timesHere],
            # igrfdf = igrf.grid(times[np.where(ind_timesHere)[0][0]],

            import igrf

            igrfdf = igrf.grid(times[ind_timesHere],
                               glat=geodesy.geodetic2geocentriclat(
                                   dfSub['gdlat'].values[ind_timesHere]),
                               glon=dfSub['gdlon'].values[ind_timesHere],
                               alt_km=dfSub['gdalt_km'].values[ind_timesHere],
                               isv=isv, itype=itype)
            
            igrfMagt = np.sqrt(igrfdf.east.values**2.+igrfdf.north.values **
                              2. + igrfdf.down.values**2.)
            
            gdlatJ, altJ, XIGRF, ZIGRF = geodesy.geoc2geod(90.-geodesy.geodetic2geocentriclat(dfSub['gdlat'].values),
                                                           geodesy.geodeticheight2geocentricR(
                                                               dfSub['gdlat'].values, dfSub['gdalt_km'].values),
                                                           -igrfdf.north.values, -igrfdf.down.values)
            EIGRF[ind_timesHere] = igrfdf.east.values
            NIGRF[ind_timesHere] = XIGRF
            UIGRF[ind_timesHere] = -ZIGRF
            igrfMag[ind_timesHere] = igrfMagt


        if return_dipoletilt:
            dptilttmp = dipole_tilt(times[ind_timesHere],epoch=datetime_to_yearfrac([apexRefTime])[0])
            dptilt[ind_timesHere] = dptilttmp

        # Increment apexRefTime by relDelta
        apexRefTime += relDelta

        nIter += 1
        if nIter >= maxIterHere:
            print("Too many iterations! Breaking ...")
            break


    if canDoMLT:

        if max_N_months_twixt_apexRefTime_and_obs == 0:

            # This can be replaced with PyAMPS code
            # /SPENCEdata/Research/Satellites/Swarm/pyAMPS/pyamps/mlt_utils.py
            # mlt = np.array([babyFunc2(a, mlonna, datime)
            #                 for mlonna, datime in zip(mlon,
            #                                           times)])
            mlt = mlon_to_mlt(mlon, times, times[0].year)

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
                    nIter += 1
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
                        mlt[tmpUseInds] = mlon_to_mlt(mlon[tmpUseInds],
                                                      times[tmpUseInds],
                                                      times[tmpUseInds][0].year)


                        nIndsConverted += nToConvert

                else:
                    mlt[ind_timesHere] = mlon_to_mlt(mlon[ind_timesHere],
                                                     times[ind_timesHere],
                                                     times[ind_timesHere][0].year)

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

    if return_dipoletilt:
        returnList.append(dptilt)
        rListNames.append('dptilt')

    if return_IGRF:

        ####################
        # OLD WAY (no loop over time)

        # isv = 0
        # itype = 1

        # glat, glon: geographic Latitude, Longitude
        # HAD TO KLUGE gridigrf12 TO GET IT TO WORK: ORIG ONLY USED ONE ALTITUDE

        # igrfdf = igrf12.gridigrf12(times,
        #                          glat=geodesy.geodetic2geocentriclat(
        #                              dfSub['gdlat'].values),
        #                          glon=dfSub['gdlon'].values,
        #                          alt_km=dfSub['gdalt_km'].values, isv=isv, itype=itype)

        # igrfMag = np.sqrt(igrfdf.east.values**2.+igrfdf.north.values **
        #                   2. + igrfdf.down.values**2.)

        # gdlatJ, altJ, XIGRF, ZIGRF = geodesy.geoc2geod(90.-geodesy.geodetic2geocentriclat(dfSub['gdlat'].values),
        #                                                geodesy.geodeticheight2geocentricR(
        #                                                    dfSub['gdlat'].values, dfSub['gdalt_km'].values),
        #                                                -igrfdf.north.values, -igrfdf.down.values)

        # output?
        # mlat,mlon,mapratio[,mlt],igrfdf.east,XIGRF(northGeodetic),-ZIGRF(upGeodetic),igrfMag

        # returnList = returnList + [igrfdf.east.values, XIGRF, -ZIGRF, igrfMag]
        # rListNames = rListNames + ['igrfE', 'igrfN', 'igrfU', 'igrfMag']

        ####################
        # NEW WAY (loop over time)

        returnList = returnList + [EIGRF, NIGRF, UIGRF, igrfMag]
        rListNames = rListNames + ['igrfE', 'igrfN', 'igrfU', 'igrfMag']

    if get_apex_basevecs:

        # if not quiet:
        #     print("Getting Apex basevectors ...")
        # # From Laundal and Richmond (2016):
        # # e1 "points eastward along contours of constant λma,"
        # # e2 "points equatorward along contours of constant φ ma (magnetic meridians)"
        # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = a.basevectors_apex(
        #     dfSub['gdlat'].values, dfSub['gdlon'].values, dfSub['gdalt_km'].values, coords='geo')

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
            returnList = returnList + [f1[0, ], f1[1, ], #f1[2, ],
                                       f2[0, ], f2[1, ], #f2[2, ],
                                       f3[0, ], f3[1, ], f3[2, ]]
            rListNames = rListNames + ['f10', 'f11', #'f12',
                                       'f20', 'f21', #'f22',
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


def get_ea_binsize_sqkm(ea):

    from hatch_python_utils.earth.geodesy import get_h2d_bin_areas
    haversine=True
    spherical_rectangle=True
    do_extra_width_calc=True
    altitude=100
    binareaopts = dict(haversine=haversine,spherical_rectangle=spherical_rectangle,do_extra_width_calc=do_extra_width_calc,altitude=altitude)
    
    EA_binarea_sqkm = get_h2d_bin_areas(ea.ea.mini, ea.ea.maxi, ea.ea.minm*15., ea.ea.maxm*15,**binareaopts)[0]
    
    #See journal__2019116__equalarea_binsize.ipynb for more info
    assert np.isclose(100187.38,EA_binarea_sqkm,rtol=1e-2),"bin size should match joinal!"

    return EA_binarea_sqkm


class EqualAreaBins(object):
    """
    See example usage in journal__20190815__lek_med_Swarm_crosstrack.ipynb

    Binsize at 0 altitude is 100187.38 sq. km
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
                        weights=None,
                        col__doNotDivideByWeight=[],
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
        candivideCols = ['j', 'je']
        if weights is not None:
            doNotDivideCols = ['deltat', 'count']

            validcols = candivideCols + doNotDivideCols
            assert all([col in validcols for col in list(data.columns)]
                       ), "Not sure how to handle some of your columns!"

            divideColInds = []
            for divideCol in candivideCols:
                candidateLoc = np.where(
                    [col == divideCol for col in list(data.columns)])[0]
                if candidateLoc.size > 0:
                    divideColInds.append(candidateLoc)
            divideColInds = np.array(divideColInds).flatten()

        if weights is None:
            def intrastatfunc(data, inds):
                return statfunc(data[inds], axis=0)
        else:
            def intrastatfunc(data, inds, weights):
                return statfunc(data[inds]*weights[inds], axis=0)

    else:
        if weights is None:
            def intrastatfunc(data, inds):
                return statfunc(data[inds, :], axis=0)
        else:
            def intrastatfunc(data, inds, weights):
                return statfunc(data[inds, :] * weights[inds, :], axis=0)

    if ea is None:
        print("Loading in equal-area thing sjølv")
        ea = EqualAreaBins().ea

    latII = pd.IntervalIndex.from_arrays(ea.mini, ea.maxi, closed='left')
    mltII = pd.IntervalIndex.from_arrays(ea.minm, ea.maxm, closed='left')

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

    statlist = np.zeros((ea.maxi.size, nCols), dtype=np.float64)*np.nan

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
            # tmpstat = statfunc(data[inds,:],axis=0)
            # tmpstat[:nCols-add_count] = intrastatfunc(data, inds)
            if weights is None:
                statlist[i, :lastdataind] = intrastatfunc(data, inds)
            else:
                statlist[i, :lastdataind] = intrastatfunc(data, inds, weights)
                # DON'T divide by weights here; typically want this to happen OUTSIDE this routine, since we divide, say, mono stats by all-time stats
                # weightsums = np.sum(weights[inds], axis=0)
                # statlist[i, divideColInds] /= weightsums.iloc[divideColInds]

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


def get_magnetic_polcap_equalarea_bin_weights(ea, apexObj,
                                              polcaplowlat=70.,
                                              ea_altitude=0,
                                              nboxes=9,
                                              mirror_SH=False,
                                              verbose=False,
                                              debug=False):

    assert nboxes == 9, "Only nboxes==9 is supported!"

    # This function tells us how to compare based on which hemisphere the ea object is set up for
    # it looks redundant to give polcaplowlat as an arg and get it back, but yes, it's intentional
    tmpcompfunc,polcaplowlat = _polcap__get_compfunc(ea,
                                                     polcaplowlat,
                                                     mirror_SH=mirror_SH)

    ea_alts = np.array([ea_altitude]*len(ea.ea.maxi))

    # weights = np.zeros((ea.maxi.size,nboxes))
    weights = np.zeros(ea.ea.maxi.size, dtype=np.float64)

    if nboxes == 9:
        ea_9sq, boxes = _get_equalarea_9square_dataframe(ea.ea)

    # breakpoint()

    for ibox, box in enumerate(boxes):
        if debug:
            print("Box #{:d}".format(ibox))

        ea_mlat, ea_mlon = apexObj.geo2apex(box[1], box[0], ea_alts)
        # incap_inds = ea_mlat >= polcaplowlat
        incap_inds = tmpcompfunc(ea_mlat, polcaplowlat)
        weights[incap_inds] += 1.

    weights = weights/nboxes

    return weights


def get_geographic_polcap_equalarea_bin_weights(ea, polcaplowlat=70.,
                                                nboxes=9,
                                                mirror_SH=False,
                                                verbose=False,
                                                debug=False):

    assert nboxes == 9, "Only nboxes==9 is supported!"

    # This function tells us how to compare based on which hemisphere the ea object is set up for
    # it looks redundant to give polcaplowlat as an arg and get it back, but yes, it's intentional
    tmpcompfunc,polcaplowlat = _polcap__get_compfunc(ea,
                                                     polcaplowlat,
                                                     mirror_SH=mirror_SH)

    # weights = np.zeros((ea.maxi.size,nboxes))
    weights = np.zeros(ea.ea.maxi.size, dtype=np.float64)

    if nboxes == 9:
        ea_9sq, boxes = _get_equalarea_9square_dataframe(ea.ea)

    # breakpoint()

    for ibox, box in enumerate(boxes):
        if debug:
            print("Box #{:d}".format(ibox))

        ea_glat, ea_glon = box[1], box[0]
        incap_inds = tmpcompfunc(ea_glat, polcaplowlat)
        weights[incap_inds] += 1.

    weights = weights/nboxes

    return weights


def get_magnetic_polcap_equalflux_bin_weights(ea, 
                                              polcaplowlat=70.,
                                              nboxes=9,
                                              mirror_SH=False,
                                              verbose=False,
                                              debug=False):

    assert nboxes == 9, "Only nboxes==9 is supported!"

    # This function tells us how to compare based on which hemisphere the ea object is set up for
    # it looks redundant to give polcaplowlat as an arg and get it back, but yes, it's intentional
    tmpcompfunc,polcaplowlat = _polcap__get_compfunc(ea,
                                                     polcaplowlat,
                                                     mirror_SH=mirror_SH)

    weights = np.zeros(ea.ea.maxi.size, dtype=np.float64)

    if nboxes == 9:
        # longitude_mode = False because we want mlts, vet du
        ea_9sq, boxes = _get_equalarea_9square_dataframe(ea.ea,
                                                         longitude_mode=False)

    for ibox, box in enumerate(boxes):
        if debug:
            print("Box #{:d}".format(ibox))

        # ea_mlat = box[1]
        incap_inds = tmpcompfunc(box[1], polcaplowlat)
        weights[incap_inds] += 1.

    weights = weights/nboxes

    return weights


def _get_equalarea_9square_dataframe(ea,
                                     longitude_mode=True):

    # Divide each into 9 boxes, count up how many
    # LON positions: LEFT, LEFT-CENTER, CENTER, CENTER-RIGHT, RIGHT [L, LC, C, CR, R]
    # LAT positions: BOTTOM, BOTTOM-MID, MID, MID-TOP, TOP          [B, BM, M, MT, T]

    # 'ea_9sq' = "Equal-area divided into nine squares"
    ea_9sq = pd.DataFrame(dict(latbm=(ea.mini+ea.centeri)/2,  
                               latm=ea.centeri,               
                               latmt=(ea.maxi+ea.centeri)/2,
                               mltlc=(ea.minm+ea.centerm)/2,
                               mltc=ea.centerm,
                               mltcr=(ea.maxm+ea.centerm)/2))

    if longitude_mode:
        ea_9sq.loc[:,'mltlc'] = ea_9sq.loc[:,'mltlc']*15.
        ea_9sq.loc[:,'mltc'] = ea_9sq.loc[:,'mltc']*15.
        ea_9sq.loc[:,'mltcr'] = ea_9sq.loc[:,'mltcr']*15.

    # Top row
    box0 = [ea_9sq['mltlc'].values, ea_9sq['latmt'].values]
    box1 = [ea_9sq['mltc'].values, ea_9sq['latmt'].values]
    box2 = [ea_9sq['mltcr'].values, ea_9sq['latmt'].values]

    # Mid row
    box3 = [ea_9sq['mltlc'].values, ea_9sq['latm'].values]
    box4 = [ea_9sq['mltc'].values, ea_9sq['latm'].values]
    box5 = [ea_9sq['mltcr'].values, ea_9sq['latm'].values]

    # Bottom row
    box6 = [ea_9sq['mltlc'].values, ea_9sq['latbm'].values]
    box7 = [ea_9sq['mltc'].values, ea_9sq['latbm'].values]
    box8 = [ea_9sq['mltcr'].values, ea_9sq['latbm'].values]

    return ea_9sq, [box0, box1, box2,
                   box3, box4, box5,
                   box6, box7, box8]


def _polcap__get_compfunc(ea,
                          polcaplowlat,
                          mirror_SH=False):

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

    return tmpcompfunc,polcaplowlat

