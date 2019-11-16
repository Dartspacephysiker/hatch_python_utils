import numpy as np


def bearing(lat1, lat2, lon1, lon2, deg=True):
    if deg:
        lat1, lat2, lon1, lon2 = map(np.deg2rad, (lat1, lat2, lon1, lon2))
    x = np.cos(lat2)*np.sin(lon2-lon1)
    y = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    return np.rad2deg(np.arctan2(x, y))


def mlt_bearing(mlats, mlts, shift_to_90deg=False):

    dlat = np.diff(mlats, append=mlats[-1]+(mlats[-1]-mlats[-2]))

    lonunwrap = np.rad2deg(np.unwrap(np.deg2rad(mlts*15.)))
    dlon = np.diff(
        lonunwrap, append=lonunwrap[-1]+(lonunwrap[-1]-lonunwrap[-2]))

    # FEIL
    # np.diff(mlts).max()
    # > 23.999415079752602
    # RIKTIG
    # np.diff(np.rad2deg(lonunwrap)/15.).max()
    # >0.27708638509115247

    lat1 = mlats
    lat2 = mlats+dlat

    lon1 = lonunwrap
    lon2 = lonunwrap+dlon

    if shift_to_90deg:
        magbearing = np.abs(bearing(lat1, lat2, lon1, lon2, deg=True))
        magbearing[magbearing > 90] -= 180
        magbearing = np.abs(magbearing)
        return magbearing
        # mag.loc[(mag['magbearingShift'] > 90),'magbearingshift'] = mag[mag['magbearingShift'] > 90]['magbearingShift']-180
        # mag.loc[mag['magbearingShift'] > 90,'magbearingShift'] -= 180
    else:
        return bearing(lat1, lat2, lon1, lon2, deg=True)
