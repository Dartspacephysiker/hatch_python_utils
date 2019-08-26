# 2019/08/26
# def Schunk_Nagy_ion_neutral_coll_freqs():
# Provide


def ion_neutral_coll_freq__H_Hplus(nn, Ti, Tn):
    """
    Applicable for Tr > 50
    """
    Tr = (Ti+Tn)/2.
    threshTr = 50               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 2.65e-10 * nn * np.sqrt(Tr) * (1-0.083*np.log10(Tr))**2


def ion_neutral_coll_freq__He_Heplus(nn, Ti, Tn):
    """
    Applicable for Tr > 50
    """
    Tr = (Ti+Tn)/2.
    threshTr = 50               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 8.73e-11 * nn * np.sqrt(Tr) * (1-0.093*np.log10(Tr))**2


def ion_neutral_coll_freq__N_Nplus(nn, Ti, Tn):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 275               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 3.83e-11 * nn * np.sqrt(Tr) * (1-0.063*np.log10(Tr))**2


def ion_neutral_coll_freq__O_Oplus(nn, Ti, Tn):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 235               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 3.67e-11 * nn * np.sqrt(Tr) * (1-0.064*np.log10(Tr))**2


def ion_neutral_coll_freq__N2_N2plus(nn, Ti, Tn):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 170               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 5.14e-11 * nn * np.sqrt(Tr) * (1-0.069*np.log10(Tr))**2


def ion_neutral_coll_freq__O2_O2plus(nn, Ti, Tn):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 800               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 2.59e-11 * nn * np.sqrt(Tr) * (1-0.073*np.log10(Tr))**2


def ion_neutral_coll_freq__O_Hplus(nn, Ti, Tn):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 300               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 6.61e-11 * nn * np.sqrt(Tr) * (1-0.047*np.log10(Tr))**2


def ion_neutral_coll_freq__H_Oplus(nn, Ti, Tn):
    """
    Applicable for Tr > threshTr
    """
    Tr = (Ti+Tn)/2.
    threshTr = 300               # K
    if Tr <= threshTr:
        print("Whoa! Below limit!")
    return 4.63e-12 * nn * np.sqrt(Tn+Ti/16.)


def ion_neutral_coll_frqe
