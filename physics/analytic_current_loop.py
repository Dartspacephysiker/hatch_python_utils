"""
Expressions taken from, or derived starting with,

Hampton, S. et al. (2020) ‘Closed-form expressions for the magnetic fields of rectangular and
circular finite-length solenoids and current loops’, AIP Advances, 10(6), p. 065320. doi: 10.1063/5.0010982.

"""
import numpy as np
mu0 = np.pi*4e-7
from scipy.special import ellipk,ellipe
from pysymmetry.constants import RE


def current_estimate_for_Dst(Dst_nT,a_RE):
    """
    Estimate the current needed to produce a particular value of Dst assuming a ring current distance of a_RE.
    This expression is obtained by starting with Equation (19) of Hampton et al (2020), then
    setting r = 1 R_E, z = 0, B0 = mu0 * I / (2 a), and solving for I.
    """
    m = 4*a_RE/(a_RE+1)**2
    Em = ellipe(m)
    Km = ellipk(m)
    return 2 * np.pi * (RE*1e3) * (Dst_nT * 1e-9) / mu0 / (Em/(a_RE-1)+Km/(a_RE+1)) #Convert RE from km to m with 1e3, Dst from nT to T with 1e-9


def effective_dipole_moment(I,r_RE,a_RE):
    """
    Obtained by starting with Equation (19) of Hampton et al (2020), then
    setting z = 0, B0 = mu0 * I / (2 a), setting the entire expression equal to
    the z component (= -theta component) of a dipole field with theta=90°, 
    B_theta(theta=90°) = mu0 M_eff/(4 * pi * r^3), and then solving for
    M_eff.

    INPUTS
    ======
    I      : Ring current magnitude [A]
             Should be negative if it corresponds to a negative Dst value, typically of order millions of amps.
    r_RE   : Distance from Earth in equatorial plane, units of Earth radii
    a_RE   : Ring current radius, units of Earth radii

    OUPUT
    =====
    M_eff  : Effective dipole moment [A/m^2]

    For reference, Earth's dipole moment (as of 2020) is M_Earth ~= 8.22 x 10^22 A/m^2
    """
    r = r_RE * RE

    m = 4*a_RE*r_RE/(a_RE+r_RE)**2
    Em = ellipe(m)
    Km = ellipk(m)

    M_eff = 2 * I * r_RE**3 * (RE*1e3)**2 * (Em/(a_RE-r_RE)+Km/(a_RE+r_RE))
    return M_eff


def _current_loop_Bfield_cylindrical_geomfactor__rcomponent(r,z,a):
    """
    From Equation (17), Hampton et al (2020)
    https://aip.scitation.org/doi/full/10.1063/5.0010982

    """

    if not hasattr(r,"__len__"):
        r = np.array([r])

    if not hasattr(z,"__len__"):
        z = np.array([r])

    zero_r = np.isclose(r,0)
    wherezero = np.where(zero_r)[0]

    if len(wherezero) == len(r):
        return r*0

    m = 4*a*r/((a+r)**2+z**2)

    Em = ellipe(m)
    Km = ellipk(m)
    
    factor1 = a * z / np.pi / r / np.sqrt((a+r)**2+z**2)
    factor2 = (a**2+r**2+z**2)/((a-r)**2+z**2) * Em - Km

    Br = factor1 * factor2

    Br[zero_r] = 0.

    return Br


def _current_loop_Bfield_cylindrical_geomfactor__zcomponent(r,z,a):
    """
    From Equation (19), Hampton et al (2020)
    https://aip.scitation.org/doi/full/10.1063/5.0010982

    INPUTS
    ========
    r      : r coordinate
    z      : z coordinate
    a      : loop radius

    """
    m = 4*a*r/((a+r)**2+z**2)

    Em = ellipe(m)
    Km = ellipk(m)
    
    factor1 = a / np.pi / np.sqrt((a+r)**2+z**2)
    factor2 = (a**2-r**2-z**2)/((a-r)**2+z**2) * Em + Km
    
    Bz = factor1 * factor2
    return Bz


def _current_loop_Bfield_cylindrical_geomfactor__zcomponent__zEq0(r,a):
    """
    From Equation (19), Hampton et al (2020)
    https://aip.scitation.org/doi/full/10.1063/5.0010982

    INPUTS
    ========
    r      : r coordinate
    z      : z coordinate
    a      : loop radius

    """
    m = 4*a*r/((a+r)**2)

    Em = ellipe(m)
    Km = ellipk(m)
    
    factor1 = a / np.pi 
    factor2 = (Em/(a-r)+Km/(a+r))
    
    Bz = factor1 * factor2

    return Bz


def current_loop_Bfield_vector(r,z,a):

    Br = _current_loop_Bfield_cylindrical_geomfactor__rcomponent(r,z,a)
    Btheta = Br * 0
    Bz = _current_loop_Bfield_cylindrical_geomfactor__zcomponent(r,z,a)

    return [Br,Btheta,Bz]
