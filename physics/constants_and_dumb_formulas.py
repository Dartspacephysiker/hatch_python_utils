########################################
# Conversions
#
# atomic mass to eV/c^2 (20190724)
# K_to_eV               (20190724)
########################################
# Formulæ
#
# gyroradius(m, vperp, B, q=sciConst.elementary_charge)
# gyroradius_temperature(m, Tperp, B, q=sciConst.elementary_charge)
# gyrofrequency_Hz(m, B, q=sciConst.elementary_charge)
# plasma_freq_Hz(m, density)
# plasma_dens_from_f_pe_Hz(f_pe)
# plasma_dens__cm3_from_f_pe_kHz(f_pe)
########################################


import scipy.constants as sciConst

########################################
# Conversions
########################################

# ........................................
# atomic mass to eV/c^2 (20190724)
# ........................................

oMass = sciConst.atomic_mass*15.999  # kg

oMass_eV_per_cSq = oMass/sciConst.eV
cSq = sciConst.speed_of_light**2.
# oMass_eV_per_cSq,oMass_eV_per_cSq*cSq

eMass_eV_per_cSq = sciConst.electron_mass/sciConst.eV
eMass_eV_per_cSq*cSq
# >> 510998.94626861025 #Akkurat som det skal!


# ........................................
# K_to_eV      (20190724): Forandre fargene til sekundær axis
# ........................................

kelvin_to_eV = sciConst.Boltzmann/sciConst.eV
# 1/kelvin_to_eV
# >> 11604.522060401008

########################################
# Formulæ
########################################


# F.Gyroradius
def gyroradius(m, vperp, B, q=sciConst.elementary_charge):
    """
    Alle enheter er SI
    m     : masse (kg)
    vperp : hastighetskomponenten som er vinkelrett på bakgrunn-feltet (m/s)
    B     : Bakgrunnfelt størrelse (T)
    q     : Ladning på partikkel (C)

    Example
    =======
    m = 2.6567e-26 #oksygen-ion-masse
    vperp = 1000               # termisk hastighet for oksygen på 2000 km-høyde
    B = 2.3e-5                 # Feltstørrelse på 2000 km-høyde
    O_gyro = gyroradius(m,vperp,B) # gyroradius i meters
    """
    return m*vperp/q/B


def gyroradius_temperature(m, Tperp, B, q=sciConst.elementary_charge):
    """
    Alle enheter er SI
    m     : masse (kg)
    Tperp : Temperaturkomponenten vinkelrett på bakgrunn-feltet (K)
    B     : Bakgrunnfelt størrelse (T)
    q     : Ladning på partikkel (C)

    Example
    =======
    m = 2.6567e-26 #oksygen-ion-masse
    vperp = 1000               # termisk hastighet for oksygen på 2000 km-høyde
    B = 2.3e-5          # Feltstørrelse på 2000 km-høyde
    O_gyro = gyroradius(m,vperp,B) # gyroradius i meters
    """
    vperp = np.sqrt(sciConst.Boltzmann*Tperp/m)
    return m*vperp/q/B


# F.Gyrofrequency
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
    O_freq = gyrofrequency_Hz(m,B) # gyrofrekvense i Hz
    """
    return q*B/m/2/np.pi


# F.PlasmaFrequency
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
    return np.sqrt(density*sciConst.elementary_charge**2./(m*sciConst.epsilon_0))/2./np.pi


def plasma_dens_from_f_pe_Hz(f_pe):
    return (2*np.pi*f_pe)*sciConst.epsilon_0*sciConst.electron_mass/(sciConst.elementary_charge**2.)

# https://library.psfc.mit.edu/catalog/online_pubs/NRL_FORMULARY_13.pdf


# F.PlasmaDens_from_PlasmaFreq
def plasma_dens__cm3_from_f_pe_kHz(f_pe):
    return (f_pe/8.98)**2       # Frequency in kHz


# F.CollisionFreqs
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
