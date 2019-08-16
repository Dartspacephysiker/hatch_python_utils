# 2019/08/16

import numpy as np


def calc_B0_phasor(omega, sigma, epsilon, mu,
                   E0=1):
    """
    Relationship between B0 and E0 (assuming E0 phase is 0 deg)

    Returns B0 phasor in T (Wb/m^2, or N/A-m).

    omega   : Angular wave frequency (rad/s, or 2 * pi * f)
    sigma   : Material conductivity  (S/m, or A/V)
    epsilon : Material permittivity  (F/m, or N/V^2)
    mu      : Material permeability  (H/m or N/A^2)

    optional
    --------
    E0      : Electric field amp     (N/C, or V/m)

    REFS
    ====
    https://en.wikipedia.org/wiki/Wave_impedance

    EXAMPLE
    =======
    import scipy.constants as sciConst
    import numpy as np
    from hatch_python_utils.physics import em_waves_formulas as hEM

    # Get B0 in nT for params at 200-km altitude (Figure 1 in Lotko & Zhang (2018))
    ########################################

    E0      = 10                      # mV/m
    f       = 1                       # Hz  
    sigmaP  = 1e-4                    # S/m (Pedersen conductivity)
    mu      = sciConst.mu_0           # H/m (4 * pi * 1e-7 )
    epsilon = 4.97e-6                 # F/m

    B0_nT = hEM.calc_B0_phasor(2*np.pi*f,sigmaP,epsilon,mu,E0=E0*1e-3)*1e9

    phase = np.angle(B0_nT,deg=True)  # -36 degrees
    amp   = np.abs(B0_nT)             # 457.9 nT

    # Get B0 in nT for params at 2000-km altitude (Figure 1 in Lotko & Zhang (2018))
    ########################################

    E0      = 10                      # mV/m
    f       = 1                       # Hz  
    sigmaP  = 1e-12                   # S/m (Pedersen conductivity)
    mu      = sciConst.mu_0           # H/m (4 * pi * 1e-7 )
    epsilon = 7.96e-9                 # F/m

    B0_nT = hEM.calc_B0_phasor(2*np.pi*f,sigmaP,epsilon,mu,E0=E0*1e-3)*1e9

    phase = np.angle(B0_nT,deg=True)  # -0.00057 degrees
    amp   = np.abs(B0_nT)             # 100.0 nT
    """

    # Slightly more illuminating form:
    # B0 = 1/c * np.sqrt(1-(1j)*sigma/epsilon/omega) * E0

    return np.sqrt((mu * (sigma+(1j)*omega*epsilon))/((1j)*omega))
