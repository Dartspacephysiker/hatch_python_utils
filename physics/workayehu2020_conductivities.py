import numpy as np
import scipy.constants as sciConst

def coll_freqs(nN2,nO2,nO,Te,Ti,Tn):
    """
    ven, vin1, vin2, vin3 = coll_freqs(nN2,nO2,nO,Te,Ti,Tn)

    This returns collision frequencies from Workayehu et al (2020) Appendix A
    0. Electron-neutral coll. freq.   , Eq (A3).
    1. Ion-neutral coll. freq. for NO+, Eq (A4)
    2. Ion-neutral coll. freq. for O2+, Eq (A5)
    3. Ion-neutral coll. freq. for O+ , Eq (A6)
    
    Input
    ====
    nN2   : neutral N2 number density (m^-3)
    nO2   : neutral O2 number density (m^-3)
    nO    : neutral O number density  (m^-3)
    
    Te    : Electron temperature (K)
    Ti    : Ion temperature      (K)
    Tn    : Neutral temperature  (K)
    
    From Appendix A of Workayehu et al (2020):
    "The electron-neutral and ion-neutral collision frequencies are calculated using the 
    formula given by Schunk and Nagy (2009) as follows, where the three most important 
    neutral species O2, N2, and O are included"
    
    Specifically, Workayehu et al (2020) use S&N2009 Equation (4.146) 
    nu_in = C_in * n_n
    for non-resonant collisions (i.e., between ion species X+ and neutral 
    species Y where X!=Y), and S&N2009 Table 4.5 for resonant collisions.

    Also: vin1 = NO+-neutral collision frequency
          vin2 = O2+-neutral collision frequency
          vin3 = O+ -neutral collision frequency

    """

    #Reduced ion-neutral temperature
    Tr = (Ti+Tn)/2

    ven = 2.33e-17 * nN2 * (1 - 1.21e-4 * Te) * Te \
	+ 1.82e-16 * nO2 * (1 + 3.6e-2 * np.sqrt(Te) ) * np.sqrt(Te) \
	+ 8.9e-17 * nO * (1 + 5.7e-4 * Te) * np.sqrt(Te)
	
    #Ion-neutral collision frequencies
    #NO+ collision frequency
    vin1 = 4.34e-16*nN2 + 4.27e-16*nO2 + 2.44e-16*nO
	
    #O2+ collision frequency
    vin2 = 4.13e-16*nN2 + 2.31e-16*nO \
	+ 2.59e-17*nO2*np.sqrt(Tr)*(1-0.073*np.log10(Tr))**2
	
    #O+ collision frequency
    vin3 = 6.82e-16*nN2 + 6.66e-16*nO2 \
	+ 3.67e-17*nO * np.sqrt(Tr)*(1-0.064*np.log10(Tr))**2
	
    return ven, vin1, vin2, vin3


def individual_coll_freqs(nN2,nO2,nO,Te,Ti,Tn):
    """
    ve, vi1, vi2, vi3 = coll_freqs(nN2,nO2,nO,Te,Ti,Tn)

    This returns collision frequencies from Workayehu et al (2020) Appendix A
    0. Electron-neutral coll. freq.   , Eq (A3).
    1. Ion-neutral coll. freq. for NO+, Eq (A4)
    2. Ion-neutral coll. freq. for O2+, Eq (A5)
    3. Ion-neutral coll. freq. for O+ , Eq (A6)
    
    Input
    ====
    nN2   : neutral N2 number density (m^-3)
    nO2   : neutral O2 number density (m^-3)
    nO    : neutral O number density  (m^-3)
    
    Te    : Electron temperature (K)
    Ti    : Ion temperature      (K)
    Tn    : Neutral temperature  (K)
    
    From Appendix A of Workayehu et al (2020):
    "The electron-neutral and ion-neutral collision frequencies are calculated using the 
    formula given by Schunk and Nagy (2009) as follows, where the three most important 
    neutral species O2, N2, and O are included"
    
    Specifically, Workayehu et al (2020) use S&N2009 Equation (4.146) 
    nu_in = C_in * n_n
    for non-resonant collisions (i.e., between ion species X+ and neutral 
    species Y where X!=Y), and S&N2009 Table 4.5 for resonant collisions.

    Also: vin1 = NO+-neutral collision frequency
          vin2 = O2+-neutral collision frequency
          vin3 = O+ -neutral collision frequency

    """

    #Reduced ion-neutral temperature
    Tr = (Ti+Tn)/2

    ve = dict(N2=2.33e-17 * nN2 * (1 - 1.21e-4 * Te) * Te,
              O2=1.82e-16 * nO2 * (1 + 3.6e-2 * np.sqrt(Te) ) * np.sqrt(Te),
              O =8.9e-17 * nO * (1 + 5.7e-4 * Te) * np.sqrt(Te))

    #NO+ collision frequency
    vi1 = dict(N2=4.34e-16*nN2,
               O2=4.27e-16*nO2,
               O=2.44e-16*nO)
	
    #O2+ collision frequency
    vi2 = dict(N2=4.13e-16*nN2,
               O2=2.59e-17*nO2*np.sqrt(Tr)*(1-0.073*np.log10(Tr))**2,
               O=2.31e-16*nO)
	
    #O+ collision frequency
    vi3 = dict(N2=6.82e-16*nN2,
               O2=6.66e-16*nO2,
	       O =3.67e-17*nO * np.sqrt(Tr)*(1-0.064*np.log10(Tr))**2)
	
    return ve, vi1, vi2, vi3

	

def sigmap_sigmah(Ne,Nis,nu_en,nu_ins,omega_ce,omega_cis,B):
    """
    Workayehu et al, Eq (A1)
    
    INPUTS
    ======
    Ne          : Electron density (m^-3)
    Nis         : List of ion densities (NO+,O2+,O+) [m^-3]
    nu_en       : Electron-neutral collision frequency [Hz]
    nu_ins      : List of ion-neutral collision frequencies (NO+,O2+,O+) [Hz]
    omega_ce    : Electron gyrofrequency [Hz]
    omega_cis   : List of ion gyrofrequencies (NO+,O2+,O+) [Hz]
    B           : Magnetic field strength [SI]

    OUTPUTS
    =======
    sigp, sigh  : Pedersen conductivity, Hall conductivity [mho/m]
    """

    assert isinstance(Nis,list) and isinstance(nu_ins,list) and (omega_cis,list),"Ion quantities (Nis, nu_ins, omega_cis) must be lists (probably just NO+, O2+, O+)"
    assert (len(Nis) == len(nu_ins)) and (len(Nis) == len(omega_cis)),"Lists are of unequal length!"

    const = sciConst.elementary_charge*Ne/B

    # Pedersen conductivity
    electron_contribP = nu_en * omega_ce / (nu_en**2+omega_ce**2)

    ion_contribP = 0
    for i in range(len(Nis)):
        Ni = Nis[i]
        nu_in = nu_ins[i]
        omega_ci = omega_cis[i]
        ion_contribP += Ni/Ne*nu_in*omega_ci/(nu_in**2+omega_ci**2)

    sigp = const * (electron_contribP+ion_contribP)

    # Hall conductivity
    electron_contribH = omega_ce**2 / (nu_en**2+omega_ce**2)
    ion_contribH = 0
    for i in range(len(Nis)):
        Ni = Nis[i]
        nu_in = nu_ins[i]
        omega_ci = omega_cis[i]
        ion_contribH += Ni/Ne*omega_ci**2/(nu_in**2+omega_ci**2)

    sigh = const * (electron_contribH-ion_contribH)  # Here we _subtract_ the contribution from ions

    return sigp,sigh

