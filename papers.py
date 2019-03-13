# 2019/03/13
# Horwitz and Zeng outflow models

import numpy as np


def hzOPlusVelocity(EFSpecDens):
    """
    Equation (3) in Horwitz and Zeng (2009)

    EFSpecDens : "Wave electric field spectral density at O^+ gyrofreqency at 1 R_E altitude in (mV/m)^2/Hz"

    Returns:
    O+ upward velocity at 3 R_E altitude, in km/s
    """
    return 10.4 + 12. * np.tanh(3. * np.sqrt(EFSpecDens)) + 36. * np.sqrt(EFSpecDens)


def hzHPlusVelocity(EFSpecDens):
    """
    Equation (4) in Horwitz and Zeng (2009)

    EFSpecDens : "Wave electric field spectral density at O^+ gyrofreqency at 1 R_E altitude in (mV/m)^2/Hz"

    Returns:
    H+ upward velocity at 3 R_E altitude, in km/s
    """
    return 24 + 50.6 * np.sqrt(EFSpecDens)


def hzOPlusDensity(EFSpecDens, sza, fe, En):
    """
    Equation (5) in Horwitz and Zeng (2009)

    EFSpecDens : "Wave electric field spectral density at O^+ gyrofreqency at 1 R_E altitude in (mV/m)^2/Hz"
    sza        : "Solar zenith angle in degrees"
    fe         : "Electron precipitation energy flux in ergs/cm^2/s"
    En         : "Characteristic energy of electron precipitation in eV"

    Returns:
    O+ density at 3 R_E altitude, in cm^-3
    """
    tanhTerm = np.tanh(9.5-sza/10.)
    fracTerm = np.pow(fe, 1.34) / (En + 20.)**(1.25)
    return (0.0265 + 4. / (np.sqrt(EFSpecDens) + 1.84) * np.exp(-0.02/(EFSpecDens + 0.006))) \
        * (0.65 + 0.50 * tanhTerm + 625. * fracTerm
           + 10. * tanhTerm * fracTerm)


def hzHPlusDensity(EFSpecDens, sza, fe, En):
    """
    Equation (6) in Horwitz and Zeng (2009)

    EFSpecDens : "Wave electric field spectral density at O^+ gyrofreqency at 1 R_E altitude in (mV/m)^2/Hz"
    sza        : "Solar zenith angle in degrees"
    fe         : "Electron precipitation energy flux in ergs/cm^2/s"
    En         : "Characteristic energy of electron precipitation in eV"

    Returns:
    O+ density at 3 R_E altitude, in cm^-3
    """
    tanhTerm = np.tanh(9.5-sza/10.)
    fracTerm = fe / (En + 25.)**(1.05)
    return (0.1 + EFSpecDens**(0.045) / (4.55 * EFSpecDens**(0.6) + 0.91)) \
        * (0.21 + 0.18 tanhTerm + 82.1 * fracTerm + 1.43 * tanhTerm * fracTerm)


def hzOPlusFlux(EFSpecDens, sza, fe, En):
    """
    Returns:
    O+ upward flux at 3 R_E altitude, in #/cm^2/s
    """
    return 1e6 * hzOPlusVelocity(EFSpecDens) * hzOPlusDensity(EFSpecDens, sza, fe, En)  # 1e6 converts km/s to cm/s


def hzHPlusFlux(EFSpecDens, sza, fe, En):
    """
    Returns:
    H+ upward flux at 3 R_E altitude, in #/cm^2/s
    """
    return 1e6 * hzHPlusVelocity(EFSpecDens) * hzHPlusDensity(EFSpecDens, sza, fe, En)  # 1e6 converts km/s to cm/s
