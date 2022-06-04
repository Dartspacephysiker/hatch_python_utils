"""
REFERENCE
=========
Strangeway, R. J. (2012) ‘The equivalence of Joule dissipation and frictional heating in the collisional ionosphere’, Journal of Geophysical Research: Space Physics. Wiley-Blackwell, 117(A2). doi: 10.1029/2011JA017302.
"""
import numpy as np
from hatch_python_utils.physics.workayehu2020_conductivities import coll_freqs, individual_coll_freqs
from hatch_python_utils.math.vectors import dotprod
AMU = 1.6605390666e-27 #in kg

neutralmass = {'O2' : 2*15.999,  # in AMU
               'N2' : 2*14.0067,
               'O'  :   15.999}
ionmass     = {'NO+':   14.0067+15.999,
               'O2+': 2*15.999,
               'O+' :   15.999}


def fric_heating_ionneutral(nN2,nO2,nO,nNOp,nO2p,nOp,Te,Ti,Tn,Ui,Un,shape=None):
    """
    Calculate frictional heating between ions and neutrals
    Equation (28) in Strangeway (2010).


    ALL QUANTITIES MUST BE GIVEN IN SI UNITS

    INPUTS
    =======
    nN2   : neutral N2 number density (m^-3)
    nO2   : neutral O2 number density (m^-3)
    nO    : neutral O  number density (m^-3)
    
    nNOp  : NO+ ion number density (m^-3)
    nO2p  : O2+ ion number density (m^-3)
    nOp   : O+  ion number density (m^-3)
    
    Te    : Electron temperature (K)
    Ti    : Ion temperature      (K)
    Tn    : Neutral temperature  (K)

    Ui    : list of 3 array-like (m/s)
          Ion velocity vector
    Un    : list of 3 array-like (m/s)
          Neutral velocity vector

    
    """

    nidict = {'NO+':nNOp,
              'O2+':nO2p,
              'O+' :nOp}

    nndict = {'O2':nO2,
              'N2':nN2,
              'O' :nO}

    # vin1 = NO+-neutral collision frequency DICTIONARY
    # vin2 = O2+-neutral collision frequency DICTIONARY
    # vin3 = O+ -neutral collision frequency DICTIONARY
    ven, vin1, vin2, vin3 = individual_coll_freqs(nN2,nO2,nO,Te,Ti,Tn)

    vindict = {'NO+':vin1,
               'O2+':vin2,
               'O+' :vin3}

    # veldiffmag = dotprod(Ui-Un,Ui-Un)
    veldiff = [Uii-Uni for Uii,Uni in zip(Ui,Un)]
    veldiffmagsq = sum([vdi**2 for vdi in veldiff])

    fricdict = {}
    for i in ['NO+','O2+','O+']:
        print("ion: ",i)
        ionfricratedict = dict()
        mi = ionmass[i]
        ni = nidict[i]
        nui = vindict[i]

        for n in ['O2','N2','O']:
            print("   neutral: ",n)
            nuin = nui[n]
            mn = neutralmass[n]
            nn = nndict[n]

            ionfricrate = mi*mn/(mi+mn) *ni* nuin * veldiffmagsq
            if shape is not None:
                ionfricrate = ionfricrate.reshape(shape)

            ionfricratedict[n] = ionfricrate

        ionfricratedict['total'] = sum([fricrate for key,fricrate in ionfricratedict.items()])

        fricdict[i] = ionfricratedict

    fricdict['total'] = sum([ionfrd['total'] for key,ionfrd in fricdict.items()])

    return fricdict
