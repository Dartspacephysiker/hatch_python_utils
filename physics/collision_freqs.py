import scipy.constants as sciConst
import numpy as np

# Originally found in Lotko_and_Zhang_2019_model_helpers.py
def ion_neutral_coll_freq_Hz(n_n, n_i, M):
    """
    Kelley (1989) Appendix B.1, page 460

    n_n   : Neutral density       (cm^-3)
    n_i   : Ion density           (cm^-3)
    M     : Mean molecular weight (atomic mass units)
    """
    return (2.6e-9)*(n_n+n_i)/np.sqrt(M)


def electron_neutral_coll_freq_Hz(n_n, T_e):
    """
    Kelley (1989) Appendix B.1, page 462

    n_n   : Neutral density       (cm^-3)
    T_e   : Electron temperature  (K)
    M     : Mean molecular weight (atomic mass units)
    """
    return (5.4e-10)*n_n*np.sqrt(T_e)


def electron_ion_coll_freq_Hz(n_e, T_e):
    """
    Kelley (1989) Appendix B.1, page 463

    n_e   : Electron density      (cm^-3)
    T_e   : Electron temperature  (K)
    """
    return (34+4.18*np.log(T_e**3. / n_e))*n_e*T_e**(-3/2.)


def gyrofrequency_Hz(m, B, q=sciConst.elementary_charge):
    """
    Alle enheter er SI
    m     : masse (kg)
    B     : Bakgrunnfelt størrelse (T)
    q     : Ladning på partikkel (C)

    Example
    =======
    m = 2.6567e-26 #oksygen-ion-masse
    B = 2.3e-5                 # Feltstørrelse på 2000 km-høyde
    O_freq = gyrofrequency(m,B) # gyroradius i meters
    """
    return q*B/(m*2.*np.pi)


def plasma_freq_Hz(m, density):
    """
    Alle enheter er SI
    m        : effektiv masse (kg)
    density  : Tetthet (m^-3)

    Example
    =======
    m = sciConst.electron_mass #9.11e-31 kg
    density=2e9 #elektrontetthet på 2000 km-høyde
    f_pe = plasma_freq_Hz(m,density) # elektron_plasma
    """
    return np.sqrt(density*(sciConst.elementary_charge**2.)/(m*sciConst.epsilon_0))/2./np.pi


def collisionfrequency__N2(nN2, Te):
    """
    Schunk and Nagy Table 4.6
    nN2     : Neutral N2 density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 2.33e-11 * nN2 * (1-1.21e-4*Te)*Te


def collisionfrequency__O2(nO2, Te):
    """
    Schunk and Nagy Table 4.6
    nO2     : Neutral O2 density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 1.82e-10 * nO2 * (1+3.6e-2*np.sqrt(Te))*np.sqrt(Te)


def collisionfrequency__O(nO, Te):
    """
    Schunk and Nagy Table 4.6
    nO      : Neutral O density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 8.9e-11 * nO * (1+5.7e-4*Te)*np.sqrt(Te)


def collisionfrequency__He(nHe, Te):
    """
    Schunk and Nagy Table 4.6
    nHe     : Neutral He density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 4.6e-10 * nHe * np.sqrt(Te)


def collisionfrequency__H(nH, Te):
    """
    Schunk and Nagy Table 4.6
    nH      : Neutral H density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 4.5e-9 * nH * (1-1.35e-4*Te)*np.sqrt(Te)


def collisionfrequency__CO(nCO, Te):
    """
    Schunk and Nagy Table 4.6
    nCO     : Neutral CO density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 2.34e-11 * nCO * (Te+165)


def collisionfrequency__CO2(nCO2, Te):
    """
    Schunk and Nagy Table 4.6
    nCO2    : Neutral CO2 density   (cm^-3)
    Te      : Electron temperature (K)
    """
    return 3.68e-8 * nCO2 * (1+4.1e-11*np.abs(4500-Te)**(2.93))
