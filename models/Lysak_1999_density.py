# 2019/07/31
# def Lysak_1999_density():

import numpy as np
import pandas as pd


def tetthet_profile(z_km,
                    dayside=True,
                    solar_max=True,
                    units='cm',
                    return_DataFrame=True):
    """
    From Lysak (1999), "Propagation of Alfvén waves through the ionosphere: Dependence on ionospheric parameters"

    z_km     : array_like, heights at which to calculate electron density
    units    : 'cm' or 'm'. You'll get density in cm^-3 and m^-3, respectively.

    nightside and ~solar_max have not been implemented!

    """
    # From Lysak (1999) Fig 1a

    doDayMax = dayside and solar_max
    doDayMin = dayside and (not solar_max)
    # doDayMin = False
    # units = 'm'

    if doDayMax:
        # Dayside max (E layer, F1 layer, F2 layer)
        plotTit = "Dayside max"
        plotCol = "orange"

        if units == 'm':
            n0 = np.array([1.5e11, 2.5e11, 2e12])  # m^-3
        else:
            n0 = np.array([1.5e5, 2.5e5, 2e6])  # cm^-3

        z0 = np.array([100, 200, 350])
        htop = np.array([300, 200, 175])
        hbot = np.array([10, 30, 75])
    elif doDayMin:
        # Dayside min (E layer, F1 layer, F2 layer)
        plotTit = "Dayside min"
        plotCol = "blue"
        if units == 'm':
            n0 = np.array([1.5e11, 2.5e11, 5e11])  # m^-3
        else:
            n0 = np.array([1.5e6, 2.5e5, 5e5])  # cm^-3

        z0 = np.array([100, 200, 350])
        htop = np.array([400, 100, 130])
        hbot = np.array([10, 30, 75])

    else:
        assert 2 < 0, "NO"

    # nLayers = []
    nTot = np.zeros(len(z_km), dtype=np.float64)
    for layerInd in range(len(z0)):
        xLayer = np.where(z_km > z0[layerInd],
                          (z_km-z0[layerInd])/htop[layerInd],
                          (z0[layerInd]-z_km)/hbot[layerInd])
        nLayer = n0[layerInd] * np.exp(1-xLayer-np.exp(-xLayer))
        nTot += nLayer
        # nLayers.append(nLayer)

    # GAMLE MÅTEN
    # layerInd = 0                    # E layer
    # xELayer = np.where((z_km >z0[layerInd]),(z_km-z0[layerInd])/htop[layerInd],(z0[layerInd]-z_km)/hbot[layerInd])
    # nELayer = n0[layerInd] * np.exp(1-xELayer-np.exp(-xELayer))

    # layerInd = 1                    # F1 layer
    # xF1Layer = np.where((z_km >z0[layerInd]),(z_km-z0[layerInd])/htop[layerInd],(z0[layerInd]-z_km)/hbot[layerInd])
    # nF1Layer = n0[layerInd] * np.exp(1-xF1Layer-np.exp(-xF1Layer))

    # layerInd = 2                    # F2 layer
    # xF2Layer = np.where((z_km >z0[layerInd]),(z_km-z0[layerInd])/htop[layerInd],(z0[layerInd]-z_km)/hbot[layerInd])
    # nF2Layer = n0[layerInd] * np.exp(1-xF2Layer-np.exp(-xF2Layer))
    # nTot = nELayer + nF1Layer + nF2Layer

    if return_DataFrame:
        return pd.DataFrame(dict(h=z_km, n_plasma=nTot))
    else:
        return nTot
