# 2019/08/26
# def Schunk_Nagy_ion_neutral_coll_freqs():
# Provide
import numpy as np
import scipy.stats as ss


belowThreshMelding = "Whoa! Below limit!"


def ion_neutral_coll_freq__H_Hplus(nn, Tn, Ti):
    """
    Applicable for Tr > 50
    """
    Tr = (Ti+Tn)/2.
    threshTr = 50               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 2.65e-10 * nn * np.sqrt(Tr) * (1-0.083*np.log10(Tr))**2


def ion_neutral_coll_freq__He_Heplus(nn, Tn, Ti):
    """
    Applicable for Tr > 50
    """
    Tr = (Ti+Tn)/2.
    threshTr = 50               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 8.73e-11 * nn * np.sqrt(Tr) * (1-0.093*np.log10(Tr))**2


def ion_neutral_coll_freq__N_Nplus(nn, Tn, Ti):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 275               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 3.83e-11 * nn * np.sqrt(Tr) * (1-0.063*np.log10(Tr))**2


def ion_neutral_coll_freq__O_Oplus(nn, Tn, Ti):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 235               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 3.67e-11 * nn * np.sqrt(Tr) * (1-0.064*np.log10(Tr))**2


def ion_neutral_coll_freq__N2_N2plus(nn, Tn, Ti):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 170               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 5.14e-11 * nn * np.sqrt(Tr) * (1-0.069*np.log10(Tr))**2


def ion_neutral_coll_freq__O2_O2plus(nn, Tn, Ti):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 800               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 2.59e-11 * nn * np.sqrt(Tr) * (1-0.073*np.log10(Tr))**2


def ion_neutral_coll_freq__O_Hplus(nn, Tn, Ti):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 300               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 6.61e-11 * nn * np.sqrt(Tr) * (1-0.047*np.log10(Tr))**2


def ion_neutral_coll_freq__H_Oplus(nn, Tn, Ti):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 300               # K
    if np.where(Tr <= threshTr)[0].size > 0:
        print(belowThreshMelding)
    return 4.63e-12 * nn * np.sqrt(Tn+Ti/16.)


def getStatDensFunc(only_neutral=True,
                    only_ion=False,
                    summum=False,
                    gm=False,
                    hm=False,
                    mean=False,
                    minimum=False,
                    maximum=False):

    if only_neutral:
        def func(array, axis=None): return array[0]
    elif only_ion:
        def func(array): return array[1]
    elif summum:
        func = np.sum
    elif gm:
        func = ss.gmean
    elif hm:
        func = ss.hmean
    elif mean:
        func = np.mean
    elif minimum:
        func = np.min
    elif maximum:
        func = np.max

    return func


def ion_neutral_coll_freqs(neutral_list, ion_list,
                           nn_list, ni_list,
                           Tn_list, Ti_list,
                           list_output=False,
                           ion_species_avg=False,
                           use_only_neutral_density=False,  # Actually True, see below!
                           use_only_ion_density=False,
                           use_sum_of_neutral_and_ion_density=False,
                           use_geometric_mean_of_neutral_and_ion_density=False,
                           use_harmonic_mean_of_neutral_and_ion_density=False,
                           use_mean_of_neutral_and_ion_density=False,
                           use_minimum_of_neutral_and_ion_density=False,
                           use_maximum_of_neutral_and_ion_density=False,
                           verbose=False):
    """
    Calculate ion-neutral collision frequencies Table 4.4 C_in coefficients fra Schunk and Nagy (2000? 2009?)
    WARNING! Formulas assume ion-neutral pairs, which is simply not the case most of the time
    """

    # if len(nn_list) > 0:
    #     checknn = nn_list[0]
    # else:
    checknn = nn_list[0]

    try:
        junk = checknn.values
        thisisPandas = True

        thisisProfile = junk.size > 1
    except:
        thisisPandas = False
        thisisProfile = False

    pickDens_kw = dict(only_neutral=use_only_neutral_density,
                       only_ion=use_only_ion_density,
                       summum=use_sum_of_neutral_and_ion_density,
                       gm=use_geometric_mean_of_neutral_and_ion_density,
                       hm=use_harmonic_mean_of_neutral_and_ion_density,
                       mean=use_mean_of_neutral_and_ion_density,
                       minimum=use_minimum_of_neutral_and_ion_density,
                       maximum=use_maximum_of_neutral_and_ion_density)

    # Sikre at vi har valgt noe, men bare ett noe
    nPickDens_kws_set = np.sum(np.array(list(pickDens_kw.values())))
    if nPickDens_kws_set == 0:
        pickDens_kw['only_neutral'] = True

    elif nPickDens_kws_set > 1:

        print("Too many keywords set!")
        for key in pickDens_kw:
            print("{:s} : ".format(key, pickDens_kw[key]))

        return None

    if verbose:
        for key in pickDens_kw:
            if pickDens_kw[key]:
                print("Using {:s} as density for ion-neutral calc".format(key))

    statFunc = getStatDensFunc(**pickDens_kw)

    collFreq = 0
    listie = []

    nu_species = dict()
    # All list entries are in this order:
    # H+, He+, C+, N+, O+, CO+, N2+, NO+, O2+, CO2+
    # H, He, N, O, CO, N2, O2, CO2

    coeffs = load_SchunkNagy_ion_neutral_coeffs()

    titleFmtStr = "{:5s}  {:5s}  {:10s}  {:10s}  {:10s}  {:10s}  {:10s}  {:10s}  {:10s}"
    printFmtStr = "{:5s}  {:5s}  {:10.4g}  {:10.4g}  {:10.4g}  {:10.4g}  {:10.4g}  {:10.4g}  {:10.4g}"
    print(titleFmtStr.format("Neutr", "Ion",
                             "Combdens", "nn", "ni", "Tn", "Ti", "nu", "Alt"))

    for ion, ni, Ti in zip(ion_list, ni_list, Ti_list):

        tmpni = ni.copy()
        fixni = np.where(tmpni <= 0)[0]
        if fixni.size > 0:
            tmpni[fixni] = 1e-10

        if ion_species_avg:
            collFreq = 0
            denom_dens = 0

        for neutral, nn, Tn in zip(neutral_list, nn_list, Tn_list):

            # pickDens(use_geometric_mean_of_neutral_and_ion_density=True,
            #          use_harmonic_mean_of_neutral_and_ion_density=True,
            #          use_minimum_of_neutral_and_ion_density=True):

            #

            # SIKRE AT ALLE VERDIENE ER > 0
            tmpnn = nn.copy()
            fixnn = np.where(tmpnn <= 0)[0]
            if fixnn.size > 0:
                tmpnn[fixnn] = 1e-10

            # HVILKEN TETTHET BRUKER VI?
            if use_only_neutral_density:
                dens = tmpnn
            elif use_only_ion_density:
                dens = tmpni
            else:
                if thisisPandas:
                    dens = statFunc(
                        np.array([tmpnn.values, tmpni.values]), axis=0)

            # HENTE UT KOEFFISIENT
            coeff = coeffs[neutral][ion]

            # REGNE UT KOLLFREK
            if callable(coeff):
                tmpCollFreq = coeff(dens, Tn, Ti)
            else:
                tmpCollFreq = coeff * 1e-10 * dens

            # HOW TO OUTPUT
            if list_output:
                listie.append((neutral, ion, tmpCollFreq))
            elif ion_species_avg:
                denom_dens += tmpnn
                collFreq += tmpnn * tmpCollFreq
            else:
                collFreq += tmpCollFreq

            # BARE NOEN DETALJER
            if verbose:
                if thisisPandas:
                    if thisisProfile:
                        tmpind = np.int64(np.argmax(tmpni))
                        # print(tmpind)
                        print(printFmtStr.format(neutral, ion,
                                                 dens[tmpind],
                                                 tmpnn.iloc[tmpind], tmpni.iloc[tmpind],
                                                 Tn.iloc[tmpind], Ti.iloc[tmpind],
                                                 tmpCollFreq[tmpind],
                                                 tmpnn.index[tmpind]))
        if ion_species_avg:
            print("Getting avg for {:s}".format(ion))
            collFreq = collFreq / denom_dens
            nu_species[ion] = collFreq

    if list_output:
        return listie
    elif ion_species_avg:
        return nu_species
    else:
        return collFreq


def load_SchunkNagy_ion_neutral_coeffs():
    """
    Load Table 4.4 C_in coefficients fra Schunk and Nagy (2000? 2009?)
    """

    ionKeys = ['H+', 'He+', 'C+', 'N+', 'O+',
               'CO+', 'N2+', 'NO+', 'O2+', 'CO2+']

    Hcoeffs = dict(zip(ionKeys,
                       [ion_neutral_coll_freq__H_Hplus, 4.71, 1.69,
                        1.45, ion_neutral_coll_freq__H_Oplus, 0.74,
                        0.74, 0.69, 0.65,
                        0.47]))
    Hecoeffs = dict(zip(ionKeys,
                        [10.6, ion_neutral_coll_freq__He_Heplus, 1.71,
                         1.49, 1.32, 0.79,
                         0.79, 0.74, 0.70,
                         0.51]))
    Ncoeffs = dict(zip(ionKeys,
                       [26.1, 11.9, 5.73,
                        ion_neutral_coll_freq__N_Nplus, 4.62, 2.95,
                        2.95, 2.79, 2.64,
                        2.00]))
    Ocoeffs = dict(zip(ionKeys,
                       [ion_neutral_coll_freq__O_Hplus, 10.1, 4.94,
                        4.42, ion_neutral_coll_freq__O_Oplus, 2.58,
                        2.58, 2.44, 2.31,
                        1.76]))
    COcoeffs = dict(zip(ionKeys,
                        [35.6, 16.9, 8.74,
                         7.90, 7.22, None,
                         4.84, 4.59, 4.37,
                         3.40]))
    N2coeffs = dict(zip(ionKeys,
                        [33.6, 16.0, 8.26,
                         7.47, 6.82, 4.24,
                         ion_neutral_coll_freq__N2_N2plus, 4.34, 4.13,
                         3.22]))
    O2coeffs = dict(zip(ionKeys,
                        [32.0, 15.3, 8.01,
                         7.25, 6.64, 4.49,
                         4.49, 4.27, ion_neutral_coll_freq__O2_O2plus,
                         3.18]))
    CO2coeffs = dict(zip(ionKeys,
                         [41.4, 20.0, 10.7,
                          9.73, 8.95, 6.18,
                          6.18, 5.89, 5.63,
                          None]))

    coeffs = dict(H=Hcoeffs,
                  He=Hecoeffs,
                  N=Ncoeffs,
                  O=Ocoeffs,
                  CO=COcoeffs,
                  N2=N2coeffs,
                  O2=O2coeffs,
                  CO2=CO2coeffs)

    return coeffs
