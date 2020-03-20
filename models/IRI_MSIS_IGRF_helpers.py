# 2019/08/22
import scipy.constants as sciConst

import apexpy
import igrf12
import numpy as np
import pandas as pd
from load_F107DB import F107DB
import pyglow
from pytt.earth import geodesy
from functools import reduce


# IRI LIMITS from iri.org:
# HEIGHT RANGE:
# Electron density: daytime: 65-2000km, nighttime: 80-2000km
# Electron and ion temperature: 60-2500km (IRI-95 option: 60-3000km)
# Ion composition: 75-2000km (DS95/DY85 option: 80-2000km)

# IRIspecies = ['nO', 'nH', 'nHE', 'nO2', 'nNO', 'nN', 'nCluster']
IRIspecies = ['nO', 'nH', 'nHE', 'nO2', 'nNO', 'nN']
IRIKeepSpecies__density = ['ne', 'ni',
                           'nO', 'nH', 'nHE', 'nO2', 'nNO', 'nN']

MSISKeepSpecies__density = ['nn', 'HE', 'O', 'N2', 'O2',
                            'AR', 'H', 'N', 'O_anomalous']


def get_IGRF_profile(picktime, glats, glons, z_km, igrf_isv, igrf_itype):
    igrf = igrf12.gridigrf12(picktime,
                             glat=glats.flatten(),
                             glon=glons.flatten(),
                             alt_km=z_km,
                             isv=igrf_isv, itype=igrf_itype)

    igrfMag = np.sqrt(igrf.east.values**2.+igrf.north.values **
                      2. + igrf.down.values**2.)

    igrfeast = igrf.east.values
    igrfnorth = igrf.north.values
    igrfdown = igrf.down.values

    igrf = pd.DataFrame(
        dict(east=igrfeast, north=igrfnorth, down=igrfdown, mag=igrfMag))
    igrf.set_index(pd.Series(z_km, name='Height'), inplace=True)

    return igrf


def get_IGRF_IRI2016_MSIS_profile(picktime, z_km, glats, glons,
                                  convert_from_apex_mlat_mlt_to_glat_glon=False,
                                  apex_mlat=None,
                                  apex_mlt=None,
                                  apex_ref_alt_km=None,
                                  get_IGRF=True,
                                  get_IRI=True,
                                  get_MSIS=True,
                                  set_neg_vals_to_zero=False,
                                  iri_use_F107=False,
                                  msis_use_F107=False,
                                  iri__const_F107_value=None,
                                  msis__const_F107_value=None,
                                  do_scale_dens_to_mneg3=False,
                                  extrapolate_model_above_max_valid_height=True,
                                  extrapolate_model_below_min_valid_height=True,
                                  extrapolate_model_below__scaleHeight_km=8,
                                  make_IRI_Te_eq_MSIS_Tn_below_and_above_valid_IRI_range=True,
                                  dataframe=False,
                                  F107=None,
                                  verbose=False):
    """
    If dataframe, return df with any of IGRF, IRI, MSIS
    Else return list with members [IGRF,IRI,MSIS] any of which could be blank
    """
    global IRIspecies, IRIKeepSpecies__density, MSISKeepSpecies__density

    iri__have_const_F107_value = iri__const_F107_value is not None
    if iri__have_const_F107_value:
        iri__const_F107_value = np.float64(iri__const_F107_value)
        if verbose:
            print("IRI  constant F10.7 value: {:.2f}".format(iri__const_F107_value))

    msis__have_const_F107_value = msis__const_F107_value is not None
    if msis__have_const_F107_value:
        msis__const_F107_value = np.float64(msis__const_F107_value)
        if verbose:
            print("MSIS constant F10.7 value: {:.2f}".format(msis__const_F107_value))

    igrf_isv = 0
    igrf_itype = 1

    if convert_from_apex_mlat_mlt_to_glat_glon:

        # print("{:s} {:6.2f} {:6.2f}".format(
        #     picktime.strftime("%Y-%m-%dT%H:%M:%S"), mlat, mlon))

        nZ = z_km.size

        glats = np.zeros(nZ, dtype=np.float64)
        glons = np.zeros(nZ, dtype=np.float64)

        a = apexpy.Apex(date=picktime, refh=apex_ref_alt_km)
        mlon = a.mlt2mlon(apex_mlt, picktime)

        if len(apex_mlat) == 1:
            for iH, height in enumerate(z_km):

                mgdlat, glon, error = a.apex2geo(apex_mlat, mlon, height)
                theta, r, _, _ = geodesy.geod2geoc(mgdlat, height, 0, 0)
                glats[iH] = 90.-theta
                glons[iH] = glon

            print("Converted {:.2f} MLat, {:.2f} MLT to (median) {:.2f} glat, {:.2f} glon".format(
                apex_mlat, apex_mlt, np.median(glats), np.median(glons)))

        else:

            for iH, (Alat, alon, height) in enumerate(zip(apex_mlat, mlon, z_km)):

                mgdlat, glon, error = a.apex2geo(Alat, alon, height)
                theta, r, _, _ = geodesy.geod2geoc(mgdlat, height, 0, 0)
                glats[iH] = 90.-theta
                glons[iH] = glon

            # print("Converted {:.2f} MLat, {:.2f} MLT to (median) {:.2f} glat, {:.2f} glon".format(
            #     apex_mlat, apex_mlt, np.median(glats), np.median(glons)))

    if make_IRI_Te_eq_MSIS_Tn_below_and_above_valid_IRI_range and (not get_MSIS) and get_IRI:
        print("Må be om MSIS for å fike IRI-elektrontemperatur!")
        return None

    nZ = z_km.size

    if (F107 is None) and iri_use_F107:
        F107 = F107DB(nDaysToAvg=27).db

        # Hente F10.7-indeks
        f107_i = np.argmin(np.abs(F107.index-picktime))
        f107, f107a = F107.loc[F107.iloc[f107_i].name, ['F107', 'F107_27DAvg']]

    # print("F10.7s : ", f107, f107a)

    igrf = None
    iri = None
    msis = None

    ########################################
    # IGRF-modell
    ########################################

    if get_IGRF:
        igrf = get_IGRF_profile(picktime, glats, glons,
                                z_km, igrf_isv, igrf_itype)

    ########################################
    # IRI-modell
    ########################################
    if get_IRI:

        loopers = zip(glats.flatten(),
                      glons.flatten(),
                      z_km.flatten())

        IRIdtype = np.dtype({'names': ('ne', 'Te', 'Ti', 'Tn_iri', 'NmF2', 'hmF2',
                                       'nO', 'nH', 'nHE', 'nO2', 'nNO', 'nN',
                                       'nCluster'),
                             'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8',
                                         'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
                                         'f8')})
        iri = np.rec.array(np.zeros(nZ, dtype=IRIdtype))

        for i, (lat, lon, alt) in enumerate(loopers):

            pt = pyglow.Point(picktime, lat, lon, alt, user_ind=iri_use_F107)

            if iri_use_F107:
                pt.f107 = f107
                pt.f107a = f107a
            elif iri__have_const_F107_value:
                pt.f107 = iri__const_F107_value
                pt.f107a = iri__const_F107_value

            # IRI-2016
            yunk = pt.run_iri()  # default year is 2016

            iri[i] = (pt.ne, pt.Te, pt.Ti, pt.Tn_iri, pt.NmF2, pt.hmF2,
                      pt.ni['O+'], pt.ni['H+'], pt.ni['HE+'], pt.ni['O2+'],
                      pt.ni['NO+'], pt.ni['N+'],
                      pt.ni['Cluster+'])

            # iri.loc[200,['nO','nH','nHE','nO2','nNO','nN']].sum()/iri.loc[200,['ne']]
        ########################################
        # Make IRI dataframe

        iri = pd.DataFrame(iri)
        iri.set_index(pd.Series(z_km, name='Height'), inplace=True)

        # This line will hurt you, don't uncomment
        # iri = iri[(iri['ne'] >= 0) & (iri['Te'] >= 0) &
        #           (iri['Ti'] >= 0) & (iri['Tn'] >= 0)]

        # Idea here is to prevent NaNs from killing your whole density profile
        if set_neg_vals_to_zero:
            iri['nO'][iri['nO'] < 0] = 0
            iri['nH'][iri['nH'] < 0] = 0
            iri['nHE'][iri['nHE'] < 0] = 0
            iri['nO2'][iri['nO2'] < 0] = 0
            iri['nNO'][iri['nNO'] < 0] = 0
            iri['nN'][iri['nN'] < 0] = 0
            iri['nCluster'][iri['nCluster'] < 0] = 0

        # iri['ni'] = iri[['nO', 'nH', 'nHE', 'nO2']].sum(axis=1)
        # iri['ni'] = iri[['nO', 'nH', 'nHE', 'nO2', 'nNO']].sum(axis=1)
        iri['ni'] = iri[['nO', 'nH', 'nHE', 'nO2', 'nNO', 'nN']].sum(axis=1)

        Omass, Hmass, HEmass, O2mass, NOmass, Nmass = 15.999, 1.00784, 4.002602, 15.999 * \
            2, 14.0067+15.999, 14.0067

        iri['Mion'] = (Omass*iri['nO']+Hmass*iri['nH']+HEmass*iri['nHE'] +
                       O2mass*iri['nO2']+NOmass*iri['nNO']+Nmass*iri['nN'])/iri['ni']

        # bad_ne = np.where(iri['ne'] < 0)[0]
        # if bad_ne.size > 0:
        #     # if useLysak_for_lower_ne:
        #     #     iri.loc[iri.iloc[bad_ne].index.values,
        #     #             'ne'] = n_e.iloc[bad_ne].values
        #     # else:
        #     iri.loc[iri.iloc[bad_ne].index.values, 'ne'] = 0

    ########################################
    # MSIS-modell
    ########################################
    if get_MSIS:

        loopers = zip(glats.flatten(),
                      glons.flatten(),
                      z_km.flatten())

        MSISdtype = np.dtype(dict(names=('Tn_msis', 'HE', 'O', 'N2', 'O2', 'AR', 'H', 'N', 'O_anomalous', 'rho'),
                                  formats=(('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))))

        msis = np.rec.array(np.zeros(nZ, dtype=MSISdtype))

        for i, (lat, lon, alt) in enumerate(loopers):

            ptMSIS = pyglow.Point(picktime, lat, lon, alt,
                                  user_ind=msis_use_F107)

            if msis_use_F107:
                ptMSIS.f107 = f107
                ptMSIS.f107a = f107a
            elif msis__have_const_F107_value:
                ptMSIS.f107 = msis__const_F107_value
                ptMSIS.f107a = msis__const_F107_value

            # MSIS-2000
            yM = ptMSIS.run_msis()

            # msis.append((ptMSIS.Tn_msis,ptMSIS.nn,ptMSIS.rho))
            msis[i] = (ptMSIS.Tn_msis, *list(ptMSIS.nn.values()), ptMSIS.rho)

        ########################################
        # Make MSIS dataframe
        msis = pd.DataFrame(msis)
        msis.set_index(pd.Series(z_km, name='Height'), inplace=True)
        msis['nn'] = msis[['HE', 'O', 'N2', 'O2',
                           'AR', 'H', 'N', 'O_anomalous']].sum(axis=1)
        msis['M'] = msis.rho*1e-3/msis['nn'] / \
            sciConst.atomic_mass  # mean molecular weight

    ########################################
    # Extend IRI model with exponential falloff up to 5000 km

    if get_IRI:
        if extrapolate_model_above_max_valid_height:
            if verbose:
                print(
                    "Extrapolating IRI n_i and n_e below and above valid IRI alt range with exponential falloff ...")

            # Electrons
            max_e_model_i, needModel_exp_cont = exp_model_continue__above(
                iri['ne'])

            iri.loc[iri.iloc[max_e_model_i:].index, 'ne'] = needModel_exp_cont

            # Ions
            max_ion_model_i, needModel_exp_cont = exp_model_continue__above(
                iri['ni'])

            iri.loc[iri.iloc[max_ion_model_i:].index,
                    'ni'] = needModel_exp_cont

            # And species
            print(
                "extrapolate above {:.1f} for ion species {:s} ...".format(
                    iri.iloc[max_ion_model_i].name,
                    ", ".join([spec.lstrip('n')+'+' for spec in IRIspecies])))
            for species in IRIspecies:
                max_ion_model_i, needModel_exp_cont = exp_model_continue__above(
                    iri[species])

                iri.loc[iri.iloc[max_ion_model_i:].index,
                        species] = needModel_exp_cont

        if extrapolate_model_below_min_valid_height:
            if verbose:
                print(
                    "Extrapolating IRI n_e below valid IRI alt range with exponential falloff ...")

            # Now lower edge
            min_e_model_i, needModel_exp_cont = exp_model_continue__below(
                iri['ne'],
                scaleHeight=extrapolate_model_below__scaleHeight_km)

            iri.loc[iri.iloc[:min_e_model_i].index, 'ne'] = needModel_exp_cont

            # Ioner
            min_ion_model_i, needModel_exp_cont = exp_model_continue__below(
                iri['ni'],
                scaleHeight=extrapolate_model_below__scaleHeight_km)

            iri.loc[iri.iloc[:min_ion_model_i].index,
                    'ni'] = needModel_exp_cont

            print(
                "extrapolate below {:.1f} for ion species {:s} ...".format(
                    iri.iloc[min_ion_model_i].name,
                    ", ".join([spec.lstrip('n')+'+' for spec in IRIspecies])))

            for species in IRIspecies:
                min_ion_model_i, needModel_exp_cont = exp_model_continue__below(
                    iri[species],
                    scaleHeight=extrapolate_model_below__scaleHeight_km)

                iri.loc[iri.iloc[:min_ion_model_i].index,
                        species] = needModel_exp_cont

    ########################################
    # Scaling
    ########################################
    if do_scale_dens_to_mneg3:
        if get_IRI:
            iri.loc[:, IRIKeepSpecies__density] *= 1e6

        if get_MSIS:
            msis.loc[:, MSISKeepSpecies__density] *= 1e6

    ########################################
    # Make Te the same as neutral temperature below ~60 km

    if make_IRI_Te_eq_MSIS_Tn_below_and_above_valid_IRI_range and get_IRI:

        if verbose:
            print("Making IRI Te eq MSIS Tn below and above valid IRI range ...")

        # iri_no_Te_inds = iri.iloc[np.where(iri.Te < 0)[0]].index.values
        iri_no_Te_inds_lower = iri.iloc[np.where(
            (iri.Te < 0) & (iri.index.values < 1000))[0]].index.values
        iri_no_Te_inds_upper = iri.iloc[np.where(
            (iri.Te < 0) & (iri.index.values >= 1000))[0]].index.values
        # iri_no_Te_inds
        # This is not a good idea for altitudes < 1000 km
        if iri_no_Te_inds_lower.size > 0:
            iri.loc[iri_no_Te_inds_lower,
                    'Te'] = msis.loc[iri_no_Te_inds_lower, 'Tn_msis'].values

        if iri_no_Te_inds_upper.size > 0:
            upper_Te = iri.iloc[iri.index.get_loc(
                np.min(iri_no_Te_inds_upper), method='nearest')-1].Te
            iri.loc[iri_no_Te_inds_upper, 'Te'] = upper_Te

    # Føye sammen forespurte profilene
    outtie = []

    if dataframe:
        if get_IGRF:
            outtie.append(igrf)
        if get_IRI:
            outtie.append(iri)
        if get_MSIS:
            outtie.append(msis)

        final = reduce(lambda left, right: pd.merge(left, right,
                                                    left_index=True,
                                                    right_index=True),
                       outtie)
        return final

    else:

        outtie.append(igrf)

        outtie.append(iri)

        outtie.append(msis)

        return outtie


def get_IGRF_IRI2016_MSIS_profiles(picktimes, z_km, mlats, mlons,
                                   get_IGRF=True,
                                   get_IRI=True,
                                   get_MSIS=True,
                                   do_average_profiles=True,
                                   do_scale_dens_to_mneg3=True,
                                   ref_alt_km=110,
                                   set_neg_vals_to_zero=False,
                                   iri_use_F107=True,
                                   msis_use_F107=False,
                                   extrapolate_model_above_max_valid_height=True,
                                   extrapolate_model_below_min_valid_height=True,
                                   extrapolate_model_below__scaleHeight_km=8,
                                   make_IRI_Te_eq_MSIS_Tn_below_and_above_valid_IRI_range=True,
                                   verbose=False):

    global IRIspecies, IRIKeepSpecies__density, MSISKeepSpecies__density

    igrf_isv = 0
    igrf_itype = 1

    F107 = F107DB(nDaysToAvg=27).db

    try:
        nTimes = len(picktimes)
    except:
        print("Have you provided a time array?")
        try:
            picktimes = [picktimes]
            nTimes = 1
        except:
            print("Didn't work to assume you only passed one time!")
            return None

    nZ = z_km.size

    igrfList = []
    iriList = []
    msisList = []
    # outtie = []

    try:
        nMLats = len(mlats)
        haveMultiMLats = True
    except:
        nMLats = 1
        mlats = [mlats]
        haveMultiMLats = False

    try:
        nMLons = len(mlons)
        haveMultiMLons = True
    except:
        nMLons = 1
        mlons = [mlons]
        haveMultiMLons = False

        maxNum = np.max([nMLats, nMLons, nTimes])

        if (nTimes != maxNum) or (nMLats != maxNum) or (nMLons != maxNum):
            print("Must provide equal numbers of stuff!")
            return None

    for picktime, mlat, mlon in zip(picktimes, mlats, mlons):

        print("{:s} {:6.2f} {:6.2f}".format(
            picktime.strftime("%Y-%m-%dT%H:%M:%S"), mlat, mlon))

        glats = np.zeros(nZ, dtype=np.float64)
        glons = np.zeros(nZ, dtype=np.float64)

        a = apexpy.Apex(date=picktime, refh=ref_alt_km)
        for iH, height in enumerate(z_km):
            # print(iH,height)

            mgdlat, glon, error = a.apex2geo(mlat, mlon, height)
            theta, r, _, _ = geodesy.geod2geoc(mgdlat, height, 0, 0)
            glats[iH] = 90.-theta
            glons[iH] = glon

        # glat, glon: geographic Latitude, Longitude
        # HAD TO KLUGE gridigrf12 TO GET IT TO WORK: ORIG ONLY USED ONE ALTITUDE
        # if get_IGRF:
        #     igrfList.append(get_IGRF_profile(picktime,glats,glons,z_km,igrf_isv,igrf_itype))

            # igrfList.append(igrf)

        if get_IGRF or get_IRI or get_MSIS:
            igrf, iri, msis = get_IGRF_IRI2016_MSIS_profile(
                picktime, z_km, glats, glons,
                get_IGRF=get_IGRF,
                get_IRI=get_IRI,
                get_MSIS=get_MSIS,
                set_neg_vals_to_zero=set_neg_vals_to_zero,
                iri_use_F107=iri_use_F107,
                msis_use_F107=msis_use_F107,
                do_scale_dens_to_mneg3=do_scale_dens_to_mneg3,
                extrapolate_model_above_max_valid_height=extrapolate_model_above_max_valid_height,
                extrapolate_model_below_min_valid_height=extrapolate_model_below_min_valid_height,
                extrapolate_model_below__scaleHeight_km=extrapolate_model_below__scaleHeight_km,
                make_IRI_Te_eq_MSIS_Tn_below_and_above_valid_IRI_range=make_IRI_Te_eq_MSIS_Tn_below_and_above_valid_IRI_range,
                F107=F107,
                verbose=verbose)

        igrfList.append(igrf)
        iriList.append(iri)
        msisList.append(msis)

    if (nTimes > 1) and do_average_profiles:
        if verbose:
            print("Averaging profiles ...")

        ########################################
        # Regne ut gjennomsnitt

        finalList = []

        # igrfInd, iriInd, msisInd = 0, 1, 2
        igrfInd = 0 if get_IGRF else -1
        iriInd = int(get_IGRF)
        msisInd = int(get_IGRF)+int(get_IRI)

        iriAvgCols = ['Te', 'Ti', 'Tn_iri']+IRIKeepSpecies__density
        msisAvgCols = ['M', 'Tn_msis'] + MSISKeepSpecies__density

        # Ta IGRF først
        if get_IGRF:

            igrfmag = igrfList[0]['mag'].copy()
            nIGRF = len(igrfList)
            nIGRFHeights = igrfmag.size
            for i in range(1, nIGRF):
                # print(i)
                if igrfList[i]['mag'].size != nIGRFHeights:
                    print("BLEKKSPRUT!")
                    break

                igrfmag += igrfList[i]['mag'].values

            igrfmag.loc[:] /= nIGRF
            finalList.append(igrfmag)

        # og så IRI
        if get_IRI:

            finalList.append(gjennomsnitt__IRI_MSIS(iriList, iriAvgCols))

        # og så MSIS
        if get_MSIS:
            finalList.append(gjennomsnitt__IRI_MSIS(
                msisList, msisAvgCols))

        final = reduce(lambda left, right: pd.merge(left, right,
                                                    left_index=True,
                                                    right_index=True),
                       finalList)
        return final
    else:
        return igrfList, iriList, msisList


def exp_model_continue__above(series, scaleHeight=None):  # ,
    # for_lower_edge=True):

    R_E = 6400

    # if for_lower_edge:
    #     min_model_i = np.max(np.where(series > 0)[0])
    # else:
    max_model_i = np.max(np.where(series > 0)[0])
    # max_model_height = series.iloc[max_model_i].name
    max_model_height = series.index[max_model_i]
    max_model = series.iloc[max_model_i]

    if scaleHeight is None:
        scaleHeight = max_model_height/4
        print("Using scaleHeight = {:.2f} km for above interp".format(
            scaleHeight))

    needModel_heights = series.iloc[max_model_i:].index.values

    needModel_r = R_E + needModel_heights
    max_model_r = max_model_height+R_E

    needModel_exp_cont = ((max_model_r)-needModel_r)/scaleHeight
    needModel_exp_cont = np.exp(needModel_exp_cont)*max_model

    return max_model_i, needModel_exp_cont


def exp_model_continue__below(series, scaleHeight=10):

    R_E = 6400

    min_model_i = np.min(np.where(series > 0)[0])

    min_model_height = series.index[min_model_i]
    min_model = series.iloc[min_model_i]

    needModel_heights = series.iloc[:min_model_i].index.values

    needModel_r = R_E + needModel_heights
    min_model_r = min_model_height+R_E

    needModel_exp_cont = np.abs(min_model_r-needModel_r)/scaleHeight
    needModel_exp_cont = np.exp((-1.)*needModel_exp_cont)*min_model

    return min_model_i, needModel_exp_cont


def gjennomsnitt__IRI_MSIS(dfList, avgCols):

    df = dfList[0].loc[:, avgCols].copy()

    n = len(dfList)
    nHeights = df.shape[0]

    for i in range(1, n):
        # print(i)
        if dfList[i][avgCols[0]].size != nHeights:
            print("BLEKKSPRUT!")
            break

        for col in avgCols:
            df.loc[:, col] += dfList[i][col].values

    for col in avgCols:
        df.loc[:, col] /= n

    return df


def median__IRI_MSIS(dfList, avgCols):

    assert 2 < 0, "UNFINISHED"

    n = len(dfList)
    nHeights = df.shape[0]
    nCols = len(avgCols)

    toBeSnitted = np.zeros((n, nHeights, nCols), dtype=np.float64)

    assert 2 < 0, "FICKIT"
    toBeSnitted = dfList[0].loc[:, avgCols].copy()

    for i in range(1, n):
        # print(i)
        if dfList[i][avgCols[0]].size != nHeights:
            print("BLEKKSPRUT!")
            break

        for col in avgCols:
            df.loc[:, col] += dfList[i][col].values

    for col in avgCols:
        df.loc[:, col] /= n

    return df
