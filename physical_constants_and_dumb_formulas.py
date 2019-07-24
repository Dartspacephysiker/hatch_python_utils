########################################
# Conversions
# atomic mass to eV/c^2 (20190724)
# K_to_eV               (20190724)
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
